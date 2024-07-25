# mlst-hash-experiment

I wanted to see if I find any hash collisions with alleles.

## Script

I made a script make-mutations.pl to run hashes

## Running the script

```bash
for i in $DBS/*.chewbbaca; do 
  perl scripts/make-mutations.pl --numcpus 8 $i 2>mutations/$(basename $i).log | gzip -c > mutations-justmd5/$(basename $i).tsv.gz & 
done
```

## Getting results

In these particular results, I edited the script to have the following output columns

1. Original allele ID
2. md5sum
3. sha1
4. sha256
5. crc32
6. md5_56 (56 byte md5sum)

```bash
cd mutations-noseq
for i in *.tsv.gz; do 
  echo $i; 
  for col in {2..7}; do 
    echo -n "  col $col -> "; 
    zcat $i | cut -f $col | sort -n --parallel 24 | uniq -dc | wc -l; 
  done; 
  echo; 
done | tee duplications.txt
```

## Making real mutations

This script aligns all the alleles of a locus and then runs a more sophisticated mutator with `seq-gen`.
The largest number of alleles in a locus is about 4 or 5 thousand.
So I am setting the number of alleles to 50k per locus.

Making the mutations:

```bash
mkdir mutations-real
\ls $DBS/Salmonella.chewbbaca/*.fasta | shuf | \
  xargs -n 1 -P 4 bash -c '
    echo $0 >&2; 
    b=$(basename $0); 
    perl scripts/make-real-mutations.pl --sequences 50000 $0 > mutations-real/$b.txt 2> mutations-real/$b.log
  '
```

Hashsum collision detection per locus:

```bash
cat SALM_15661.fasta.txt | \
  grep . | \
  perl -MDigest::MD5=md5_hex -lane '
    print md5_hex($_);
  ' | sort | uniq -dc
# prints nothing if no collisions
# but prints counts of collisions if so
```

Bulk detection with one hashing algo

```bash
for i in *.fasta.txt; do 
  echo $i; 
  cat $i | grep . | sort | uniq | \
    perl -MDigest::MD5=md5_hex -lane '
      print md5_hex($_);
    ' | sort | uniq -dc | sort -nr | head; 
done;
```

Putting all the hashsums into one file

```bash
for i in *.fasta.txt; do 
  echo $i >&2;
  cat $i | grep . | sort | uniq | \
    perl -MMath::BigInt -MDigest::MD5=md5_hex -MDigest::SHA=sha1_hex,sha256_hex -MArchive::Zip=computeCRC32 -lane '
      BEGIN{
        my $zeroBits = Math::BigInt->new(0);
        sub reduceMd5 {
            my ($md5, $max_bits_in_result) = @_;

            my $p = Math::BigInt->new(1) << $max_bits_in_result;
            $p -= 1;
            my $rest = Math::BigInt->new("0x$md5");
            my $result = $zeroBits->copy();

            while ($rest != 0) {
                $result = $result ^ ($rest & $p);
                $rest = $rest >> $max_bits_in_result;
            }

            return "$result";
        }
      }
      # Capture the DNA sequence in $d to make it easier to type
      $d=$F[0];
      next if($seen{$d}++); 
      print join("\t", "'$i'",$d, computeCRC32($d), md5_hex($d), reduceMd5(md5_hex($d), 56), sha1_hex($d), sha256_hex($d));
    ' 
  done | gzip -c > hashsums.tsv.gz
```

Putting hashsums and sequences together into the same tsv file

```bash
zcat hashsums.tsv.gz | \
  perl -lane '
    # In this example, the hashsum is field number 2. Change accordingly.
    my($file, $seq, $hashsum) = @F[0,1,2];
    next if($seen{$file}{$seq}++); 
    push(@{ $dup{$file}{$hashsum}}, $seq); 
    END{
      print STDERR "Found hashsums from ".scalar(keys(%dup))." files.";
      %seen=(); # clear some memory
      print STDERR "Pringing any duplicates";
      while(my($file, $hashes)=each(%dup)){
        while(my($hash, $seqs)=each(%$hashes)){
          next if(@$seqs < 2); 
          print join("\t", $file, $hash, @$seqs);
        } 
      }
    }
  ' | gzip -c > dups.tsv.gz
```

Check on a given alignment of sequences of the same checksum

```bash
zcat dups.tsv.gz | \
  tail -n +6 | \
  perl -lane '
    ($file, $hashsum) = splice(@F, 0, 2); 
    for(my $i=0;$i<@F;$i++){
      $j=$i+1; 
      print ">$j\n$F[$i]";
    } 
    last;
  ' | mafft - | goalign reformat clustal | less
```
