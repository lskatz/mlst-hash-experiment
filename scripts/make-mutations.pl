#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;

use Digest::MD5 qw/md5_base64 md5_hex/;
use Digest::SHA qw/sha1_base64 sha256_base64/;
use String::CRC32 qw/crc32/;
use Math::BigInt;

use threads;
use Thread::Queue;
use threads::shared;

use version 0.77;
our $VERSION = '0.1.1';

my @NT = qw(A C G T);
my $printLock :shared;

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help numcpus=i)) or die $!;
  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus} ||= 1;

  for my $db(@ARGV){
    mutateDatabase($db, $settings);
  }

  return 0;
}

sub mutateDatabase{
  my($db, $settings) = @_;

  my $seqQ   = Thread::Queue->new;
  my $printQ = Thread::Queue->new;

  # Initialize threads
  my @thr;
  for my $i(0..$$settings{numcpus}-1){
    $thr[$i] = threads->create(\&hashsumWorker, $seqQ, $printQ, $settings);
  }
  my $printerThread = threads->create(\&printWorker, $printQ, $settings);

  local $/ = "\n>";
  for my $file(glob("$db/*.fa $db/*.fasta")){
    logmsg "Reading from $file ...";
    open(my $fh, $file) or die "ERROR: could not open $file: $!";
    while(my $entry=<$fh>){
      chomp($entry);

      my($defline, $sequence) = split(/\n/,$entry,2);
      $sequence =~ s/\s+//g;

      my $variants = mutateSeq($sequence, 11, $settings);

      $seqQ->enqueue([$defline, $sequence, @$variants]);
    }
    close $fh;
  }

  # Terminate seqQ theads and printer thread
  logmsg "Waiting for threads to finish";
  $seqQ->enqueue(undef) for(@thr);
  $_->join for(@thr);

  logmsg "Waiting for printer thread to finish";
  $printQ->enqueue(undef);
  $printerThread->join;
}

sub hashsumWorker {
  my($seqQ, $printQ, $settings) = @_;

  # Keep track of seen sequences
  # Even though we'll keep track of seen sequences in printWorker(),
  # This will help us save time by not recalculating hashsums 
  # within the same hashsumWorker().
  my %seen;

  my @buffer;

  while(defined(my $seqArr = $seqQ->dequeue)){
    my $defline = shift(@$seqArr);
    for my $seq(@$seqArr){
      next if($seen{$seq}++);
      
      my %hashsum = (
        #md5       => md5_base64($seq),
        #sha1      => sha1_base64($seq),
        #sha256    => sha256_base64($seq),
        #crc32     => crc32($seq),
        md5_56    => reduceMd5(md5_hex($seq), 56),
        defline   => $defline,
        seq       => $seq,
      );

      push(@buffer, \%hashsum);
    }
  
    if(@buffer > 100000){
      $printQ->enqueue(\@buffer);
      @buffer = ();
    }
  }
  $printQ->enqueue(\@buffer) if(@buffer);
}

sub printWorker{
  my($printQ, $settings) = @_;

  my %seen;

  my @fields = qw(defline md5_56);
  #my @fields = qw(defline seq md5 sha1 sha256 crc32 md5_56);

  while(defined(my $hashsumArr = $printQ->dequeue)){
    lock($printLock);
    for my $hashsum(@$hashsumArr){
      next if($seen{$$hashsum{seq}}++);
      print join("\t", map{$$hashsum{$_}} @fields) . "\n";
    }
  }

}

sub mutateSeq{
  my($sequence, $reps, $settings) = @_;

  my @variant;

  my $seqLength = length($sequence);

  # Do this 11 times, but one of them will have no mutation
  for my $i(1..$reps){ 
    my $seqCopy=$sequence; # make a copy of the sequence

    # introduce a single mutation
    substr($seqCopy,rand($seqLength),1) = $NT[rand(4)];

    push(@variant, $seqCopy);

    
  }

  return \@variant;

}

# Reduces the given md5 hash to a 56 bit hash or whatever other max bits.
sub reduceMd5 {
    my ($md5, $max_bits_in_result) = @_;

    my $p = Math::BigInt->new(1) << $max_bits_in_result;
    $p -= 1;
    my $rest = Math::BigInt->new("0x$md5");
    my $result = Math::BigInt->new(0);

    while ($rest != 0) {
        $result = $result ^ ($rest & $p);
        $rest = $rest >> $max_bits_in_result;
    }

    return "$result";
}


sub usage{
  print "$0: prints hashsums for a given folder of fasta files
  Usage: $0 [options] folder1 [folder2...]
  --numcpus  1
  --help        This useful help menu
  \n";
  exit 0;
}
