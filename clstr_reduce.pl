#!/usr/bin/perl


$file90 = shift;
$segs = shift;
$reduce_rate = shift;

@clstr_nos = ();
@segs = split(/,/, $segs);
@segs_no = ();
@no2seg_idx = ();
for ($i=0; $i<@segs; $i++) {
  $seg = $segs[$i];
  if ($seg =~ /-/) { ($b, $e) = split(/-/, $seg); }
  else             {  $b= $e  = $seg; }
  for($j=$b; $j<=$e; $j++) { $no2seg_idx[$j]=$i; }
  $segs_no[$i] = 0;
}


open(TMP, $file90) || die "Can not open file";
$readin = 0;
$this_no = 0;
$this_clstr = "";
$ll;
while($ll=<TMP>) {
  if ($ll =~ /^>/ ) {
    if ($readin) { 
      $this_seg = $no2seg_idx[$this_no];
      if (($segs_no[$this_seg] % $reduce_rate) == 0) {
        print $this_clstr;
      }
      $segs_no[$this_seg]++;
    }
    $this_no=0;
    $this_clstr = $ll;
  }
  else {
    $this_clstr .= $ll;
    $readin = 1; $this_no++;
  }
}
close(TMP);

    if ($readin) { 
      $this_seg = $no2seg_idx[$this_no];
      if (($segs_no[$this_seg] % $reduce_rate) == 0) {
        print $this_clstr;
      }
      $segs_no[$this_seg]++;
    }
