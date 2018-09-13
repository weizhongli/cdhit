#!/usr/bin/perl
#
use Getopt::Std;
getopts("s:S:o:f:j:",\%opts);

die usage() unless ($opts{s} or $opts{S});

my $output            = $opts{o}; 
   $output            = "pooled" unless ($output);
my $sample_in         = $opts{s};
my $sample_command_in = $opts{S}; #### ',' delimited samples, ':' delimited entries, e.g. sample1:R1.fq:R2.fq;sample2:R1.fq:R2.fq   or sample1;sample2;sample3
my $file_list         = $opts{f};
my @file_list = qw/seq.99f seq.99f.2 seq.99f-all.clstr chimeric-small-clusters-list.txt/;
   @file_list = split(/,/, $file_list) if ($file_list); 

my $job               = $opts{j};
   $job = "otu" unless ($job);


my ($i, $j, $k, $cmd);
$cmd = `mkdir $output` unless (-e $output);

foreach $i (@file_list) {
  if (-e "$output/$i") {
    die "output dir $output & file $output/$i already exist, please remove all files from $output and re-run\n";
  }
}

######## parse NGS_samples
my @NGS_samples = ();
if (defined($sample_in)) {
  open(TMP, $sample_in) || die "can not open $sample_in";
  while($ll=<TMP>){
    next if ($ll =~ /^#/);
    next unless ($ll =~ /^\w/); chop($ll);
    my ($id, @data) = split(/\s+/,$ll);
    push(@NGS_samples, $id);
  }
  close(TMP);
}
elsif (defined($sample_command_in)) {
  my @lls = split(/,/, $sample_command_in);
  foreach $ll (@lls) {
    my ($id, @data) = split(/:/, $ll);
    push(@NGS_samples, $id);
  }
}
else {
  die "no input samples";
}

foreach $i (@file_list) {
  my $target = "$output/$i";
  foreach $j (@NGS_samples) {
    my $source = "$j/$job/$i";
    if (-e $source) {
      print STDERR "cat $source >> $target\n";
      $cmd = `cat $source >> $target`;
    }
    else {
      print STDERR "Warning, $source missing\n";
    }
  }
}

sub usage {
<<EOD;
    $0 -s sample_file -o output_dir -j job -f list_files
    -s sample data file, required unless -S is present
       File format example
#Sample data file example, TAB or space delimited for following lines
Sample_ID1 sample_data_0 sample_data_1
Sample_ID2 sample_data_0 sample_data_1
Sample_ID3 sample_data_0 sample_data_1

    -S sample data from command line, required unless -s is present
           format: Sample_ID1:sample_data_0:sample_data_0:sample_data_1,Sample_ID2:sample_data_0:sample_data_1

    -j job name
    -f list of files (delimited by ,) to pool (cat) 
EOD
}

