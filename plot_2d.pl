#!/usr/bin/perl

use Image::Magick;

$file90 = shift;
$segs = shift;
$out_file  = shift;

my $image_xx = 800;
my $image_yy = 600;
my $border   = 20;
my $border2  = $border/2;
my $lab_margin = 30;
my $scale    = 200/10;
my @seg_colors = qw//;
my @seg_sizes = qw/1 2 5 6 8 10 12 14 16 18 20/;
my $p = {
  magick     => Image::Magick->new(),
  border     => 2,
  grid_color => "#eeeeee",
  gif_file   => $out_file,
  height     => $image_yy+$border+$lab_margin,
  width      => $image_xx+$border,
};
$p->{magick}->Set(size => "$p->{width}x$p->{height}");
$p->{magick}->Read('xc:white');

my @clstr_nos = ();
open(TMP, $file90) || die "Can not open file";
$readin = 0;
my $this_no = 0;
while(my $ll=<TMP>) {
  if ($ll =~ /^>/ ) {
    if ($readin) { $clstr_nos[$this_no]++; }
    $this_no=0;
  }
  else {
    $readin = 1; $this_no++;
  }
}
close(TMP);
if ($readin) { $clstr_nos[$this_no]++; }

@segs = split(/,/, $segs);
@segs_no = ();
for ($i=0; $i<@segs; $i++) {
  $seg = $segs[$i];

  if ($seg =~ /-/) { ($b, $e) = split(/-/, $seg); }
  else             {  $b= $e  = $seg; }

  $no1 = 0; for($j=$b; $j<=$e; $j++) { $no1+= $clstr_nos[$j]; }
  print "$seg\t$no1\n";
  $no2 = int($no1/$scale); $no2++ if ($no1 % $scale);
  $segs_no[$i] = $no2;
}

my @covered = ();

my $xmin = $border2;
my $xmax = $image_xx+$border2;
my $ymin = $border2;
my $ymax = $image_yy+$border2;

$p->{magick}->Draw(primitive => 'line', points => "$xmin,$ymin,$xmax,$ymin");
$p->{magick}->Draw(primitive => 'line', points => "$xmin,$ymin,$xmin,$ymax");
$p->{magick}->Draw(primitive => 'line', points => "$xmax,$ymin,$xmax,$ymax");
$p->{magick}->Draw(primitive => 'line', points => "$xmin,$ymax,$xmax,$ymax");

my $no_segs = $#segs+1;
for ($i=$#segs; $i>=0; $i--) {
  $size = $seg_sizes[$i]-1;
  $c    = "#000000";
  $no2 = $segs_no[$i];
  for ($j=0; $j<$no2; $j++) {

    $x1 = rand($image_xx)+$border2; $x2 = $x1 + $size;
    $y1 = rand($image_yy)+$border2; $y2 = $y1 + $size;

    if ($x2 > $xmax) {$x2 = $xmax; $x1 = $x2-$size;}
    if ($y2 > $ymax) {$y2 = $ymax; $y1 = $y2-$size;}
    
    if ($size) {
      $p->{magick}->Draw(primitive => 'Rectangle',
                  points    => "$x1,$y1,$x2,$y2", fill      => $c);
    }
    else {
      $p->{magick}->Draw(primitive => 'point',
                  points    => "$x1,$y1", fill      => $c);
    }
  }
  $x1 = $border2 + int($image_xx/$no_segs)*$i;      $x2 = $x1 + $size;
  $y1 = $ymax + int($lab_margin/2) - int($size/2);  $y2 = $y1 + $size;

  if ($size) {
    $p->{magick}->Draw(primitive => 'Rectangle',
                points    => "$x1,$y1,$x2,$y2", fill      => $c);
  }
  else {
    $p->{magick}->Draw(primitive => 'point',
                points    => "$x1,$y1", fill      => $c);
  }
}

$p->{magick}->Write($p->{gif_file});







#build color index
my @colors = ();
for ($i=0; $i<$steps; $i++) {
  $j = $i/$steps;
  my ($r, $g, $b) = get_gradient_color_255_old($j);
  my $c1 = uc(sprintf("#%2x%2x%2x",$r,$g,$b));
     $c1 =~ s/ /0/g;
  $colors[$i] = $c1;
}





# return color
sub get_gradient_color_255_old {
  my $ratio = shift;
  my ($r, $g, $b);

  if ($ratio >= 0 and $ratio < 0.2 ) {
    $r = 255;
    $g = int( 1275*$ratio );
    $b = 0;
  }
  elsif ($ratio >= 0.2 and $ratio < 0.4 ) {
    $r = int ( 255 - 1275*($ratio-0.2) );
    $g = 255;
    $b = 0;
  }
  elsif ($ratio >= 0.4 and $ratio < 0.6 ) {
    $r = 0;
    $g = 255;
    $b = int ( 1275 * ($ratio-0.4) );
  }
  elsif ($ratio >= 0.6 and $ratio < 0.8 ) {
    $r = 0;
    $g = int ( 255 - 1275*($ratio-0.6) );
    $b = 255;
  }
  elsif ($ratio >= 0.8 and $ratio <= 1.0 ) {
    $r = int ( 1275 *($ratio-0.8) );
    $g = 0;
    $b = 255;
  }
  return ($r,$g,$b);
}

