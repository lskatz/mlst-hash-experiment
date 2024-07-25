#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename dirname/;
use File::Copy qw/cp mv/;
use File::Temp qw/tempdir/;
use File::Which qw/which/;

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help numcpus=i tempdir=s sequences=i)) or die $!;
  usage() if($$settings{help} || !@ARGV);
  $$settings{tempdir} //= tempdir("make-models-XXXXXX", TMPDIR => 1, CLEANUP => 1);
  $$settings{numcpus} ||= 1;
  $$settings{sequences} ||= 1;

  for my $exe(qw(goalign seq-gen modeltest-ng mafft)){
    my $path = which($exe)
      or die "ERROR: could not find $exe in your PATH";
  }

  for my $fasta(@ARGV){
    logmsg $fasta;
    my $files = makeModel($fasta, $settings);
    my $m = parseModel($$files{model});
    my $s = alignmentStats($$files{aln}, $settings);

    seqGen($fasta, $files, $m, $s, $settings);
  }

  return 0;
}

sub seqGen{
    my($fasta, $files, $m, $s, $settings) = @_;

    # Create the modified phylip/newick file as described here:
    # https://snoweye.github.io/phyclust/document/Seq-Gen.v.1.3.2/Seq-Gen.Manual.html
    # https://github.com/rambaut/Seq-Gen/blob/master/examples/seqs%2Btrees.phy
    my $treefile = "$$settings{tempdir}/".basename($fasta).".phylip.tree";
    open(my $fh, ">",$treefile) or die "ERROR: could not write to $treefile: $!";
    open(my $fh2, $$files{phy}) or die "ERROR: could not read $$files{phy}: $!";
    
    # Update the phylip header to include the correct sequence length for seq-gen
    my $phylipHeader = <$fh2>;
    $phylipHeader =~ s/^\s+|\s+$//g; # whitespace trim
    my($numTaxa, $numSites) = split(/\s+/, $phylipHeader);
    my $firstEntry = <$fh2>;
    my ($firstId, $firstSeq) = (substr($firstEntry, 0, 10), substr($firstEntry, 10));
    # Remove whitespace and dashes for the first sequence to get a good seq length
    $firstSeq =~ s/[\s+\-+]//g;
    my $seqLength = length($firstSeq);
    # Now print the header with the "correct" information: num taxa and sequence length
    print $fh "    $numTaxa $seqLength\n";

    # Print the rest of the entries with no gaps
    print $fh sprintf("%-10s%s\n", $firstId, $firstSeq);
    while(<$fh2>){
        # whitespace trim
        s/^\s+|\s+$//g;
        my($id, $seq) = (substr($_, 0, 10), substr($_, 10));
        $seq =~ s/[\s+\-+]//g;
        print $fh sprintf("%-10s%s\n", $id, $seq);
    }
    close $fh2;
    print $fh "1\n"; # just the first tree
    open(my $fh3, $$files{tree}) or die "ERROR: could not read $$files{tree}: $!";
    while(<$fh3>){
      print $fh $_;
    }
    close $fh3;
    close $fh;

    # Parameters for seq-gen
    my $freqs = join(" ", @{$$m{frequencies}});
    # note: rates cannot be applied at the same time as a general model
    my $rates = join(" ", @{$$m{"subst. rates"}});
    my $gamma = $$m{"gamma shape"};
    $gamma = 1 if($gamma !~ /\d/);

    # If we do have rates, then make the model just GTR instead of "best model"
    if($rates){
      $$m{model} = "GTR";
    }
    
    my $seqgenLog = "$treefile.seq-gen.log";
    my $numSequences = $$settings{sequences};
    my $cmd = "seq-gen -m $$m{model} -r $rates -f $freqs "
            . "-s $gamma -n $numSequences -p 1 -k 1 "
            . "$treefile 2>> $seqgenLog "
            . "| goalign reformat phylip --phylip --one-line | "
            . "tail -n +2 | cut -c 11- | grep -v '[0-9]' | grep -m $numSequences . | "
            . "sed 's/\\s//g' ";
    # Record the command before running the command
    open(my $logFh, ">", $seqgenLog) or die "ERROR: could not write to $seqgenLog: $!";
    print $logFh "Command: $cmd\n";
    close $logFh;
    logmsg "CMD $cmd";
    system($cmd);
    my $exit_code = $? >> 8;
    die "ERROR with seq-gen: $! (exit code: $exit_code)\n" . `cat $seqgenLog` if $exit_code;
}

sub alignmentStats{
    my ($aln, $settings) = @_;

    my %stat;

    # Get a reference sequence
    $stat{ref} = do{
        local $/ = "\n>";
        open(my $fh, $aln) or die "ERROR: could not open $aln: $!";
        my $entry = <$fh>;
        chomp($entry);
        
        my($defline, $seq) = split(/\n/, $entry, 2);
        $seq =~ s/[\s+\-]//g;
        $seq
    };
    
    my @line = `goalign stats -i $aln`;
    chomp(@line);

    my @header;
    for my $line(@line){
        my @F = split(/\t/, $line);
        if(@F == 2){
            $stat{$F[0]} = $F[1];
        }
        elsif(@F == 3){
            if(!@header){
                @header = split(/\t/, $line);
                next;
            }
            for(my $i=1; $i<@header;$i++){
                $stat{$header[$i]}{$F[0]} = $F[$i];
            }
        }
    }
    return \%stat;
}

sub makeModel{
    my($fastaOrig, $settings) = @_;
    
    # Put this fasta into a temp directory because there are lots of files
    my $name = basename($fastaOrig, qw(.fasta .fa .fna .fas));
    my $fasta = "$$settings{tempdir}/$name/$name.fasta";
    mkdir dirname($fasta);
    cp($fastaOrig, $fasta) or die "ERROR: could not cp $fastaOrig to $fasta: $!";

    my $alnFasta = "$fasta.aln";
    my $alnPhy   = "$fasta.phy";
    my $tree     = "$alnPhy.tree";
    my $model    = "$alnPhy.out";

    # Make an alignment from the fasta
    if(!-e $alnFasta){
        system("mafft --auto --maxiterate 100 --thread $$settings{numcpus} $fasta > $alnFasta.tmp 2> $alnFasta.log");
        die "ERROR with mafft: $!" if $?;

        # Need to rename deflines to have at most 10 characters
        open(my $fh, "$alnFasta.tmp") or die "ERROR: could not open $alnFasta.tmp: $!";
        open(my $out, ">", "$alnFasta.renamed") or die "ERROR: could not write to $alnFasta.renamed: $!";
        my $seqCounter = 0;
        while(my $line = <$fh>){
            if($line =~ /^>/){
                $seqCounter++;
                $line = ">seq$seqCounter\n";
            }
            # Sequences need to be uppercase for seq-gen
            else {
                $line = uc($line);
            }
            print $out $line;
        }
        close $out;
        close $fh;    

        mv("$alnFasta.renamed", "$alnFasta")
            or die "ERROR: could not move $alnFasta.renamed to $alnFasta: $!";
    }
    
    # convert to phylip
    if(!-e $alnPhy){
        system("goalign reformat phylip -i $alnFasta --output-strict --one-line > $alnPhy.tmp 2> $alnPhy.log");
        die "ERROR with goalign: $!" if $?;
        mv("$alnPhy.tmp", "$alnPhy")
          or die "ERROR: could not move $alnPhy.tmp to $alnPhy: $!";
    }

    # Find the evolutionary model parameters and get a basic tree
    if(!-e $model){
        system("modeltest-ng -i $alnPhy -d nt -p $$settings{numcpus} -T raxml 1> $model.stdout 2> $model.stderr");
         die "ERROR: modeltest-ng failed on $alnPhy: $! (exit code: $?)" if $?;
    }

    return {
        aln => "$alnFasta",
        phy => "$alnPhy",
        tree => "$tree",
        model => "$model",
    };
}

sub parseModel{
    my($model) = @_;

    # It's easier to just have the model in memory and parse it
    # and so let's just load it into memory
    open(my $fh, $model) or die "ERROR: could not open $model: $!";
    my $content = eval{
        local $/ = undef;
        <$fh>
    };
    close $fh;

    # Find which is the best model
    # e.g., HKY+G4, F81, etc
    my %bestModel; # count votes for best model
    while($content =~ /Best model.+\n\-+\n.+?\s+(\S+)/g){
        $bestModel{$1}++;
    }
    my $bestModel = (sort {$bestModel{$b} <=> $bestModel{$a}} keys(%bestModel))[0];

    # Find the model parameters, given the best model
    my %model;
    my @line = split(/\n/, $content);
    for(my $i=0; $i<@line; $i++){
        if($line[$i] =~ /Best model/){
            # Continue to the next line and then 
            # burn the dashes line
            $i += 2;

            # Read this section until the next dash line
            while($line[$i] !~ /^\-+/){
                my($key, $value) = split(/:/, $line[$i], 2);
                $key = lc($key);
                $value =~ s/^\s+|\s+$//g;

                $model{$key} = $value;
                $i++;
            }

            if($model{model} ne $bestModel){
                logmsg "Model was read as not the best model ($bestModel) but instead as $model{model}. Refreshing.";
                %model = ();
            }
            else{
                logmsg "Model was found: $bestModel";
                last;
            }
            
        }
    }

    # If we don't have a model yet, give up and return
    if(!keys(%model)){
        logmsg "Count not find the best model";
        return \%model;
    }

    $model{frequencies} = [split(/\s+/, $model{frequencies})];
    $model{"subst. rates"} = [split(/\s+/, $model{"subst. rates"})];

    return \%model;
}

sub usage{
    my($settings) = @_;
    print "Usage: $0 [options] file.fasta [file2.fasta ...]
    --numcpus   1
    --sequences 1  How many sequences to generate
    --tempdir  
    --help
    \n";
    exit 0;
}