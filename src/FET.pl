#!/usr/bin/perl

use Storable;
use strict;
use Text::NSP::Measures::2D::Fisher::right;

my $clstr_file = shift;
my $anno_file = shift;
my $store_file = shift;

my @cls_list = ();
my @fun_list = ();
my $cur_cls = "";
my %cls2rep = (); 
my @cur_anno = ();


open(TMP, $clstr_file) || die;
while(my $ll = <TMP>) { # read .clstr files
    if ($ll =~ /^>/) { # the begin of a cluster
        $cur_cls = $ll;
        $cur_cls =~ s/^>(.*?)\s$/$1/;
    #    print "$cur_cls|\n";
    }
    else{
        chop($ll);
        if ($ll =~ /^(\d+)\s+(\d+)(aa|nt),\s+>(.+)\.\.\./) {
            my @tmp = split(/\|\|/,$4);
            if ($#tmp == 0){
                @cur_anno = ();
            }
            else{
                @cur_anno = split(/,/, pop(@tmp));
            }
#            print $cur_cls.$cur_anno[0]."|\n";
            push(@cls_list, $cur_cls);
            push(@fun_list, [@cur_anno]);
            if ($ll =~ /^(\d+)\s+(\d+)(aa|nt),\s+>(.+)\.\.\.(.*)\*$/){
     #           print "$4\n";
                $cls2rep{$cur_cls} = $4;
#                print "$cur_cls\t$4\n";
            }
        }
    }
}

#print join("\n", @cls_list[0..10]);
@cls_list = map {$cls2rep{$_}} @cls_list;
#print join("\n", @cls_list[0..10]);
#print "\n";
#foreach my $i (0..10){
#    print join("\t",@{$fun_list[$i]});
#    print "\n";
#} 
#print join("\n", @fun_list[0..10]);
#exit(1);
my %cls_size = ();
my %cls_anno = ();
my %anno_size = ();
my $M = $#fun_list+1;
#print $#fun_list."\t".$M."\n";
#print $#cls_list."\t".$M."\n";
foreach my $i (0..$#fun_list){
    $cls_size{$cls_list[$i]}++;
    if ($#{$fun_list[$i]} >= 0) { # have annotation
        foreach my $anno (@{$fun_list[$i]}){
#            print "$i\t$cls_list[$i]\t$anno\n";
            $anno_size{$anno}++;
            $cls_anno{$cls_list[$i]}{$anno}++;
        }
    }
}

#while (my ($a,$b) = each %anno_size){
#    print "$a\t$b\n";
#}

#print "COG0171\t".$anno_

my %resu = ();
while(my ($cls, $cls_s) = each %cls_size){
    my @tmp = ();
#    $resu{$cls} = [];
    while (my ($anno,$anno_s) = each %{$cls_anno{$cls}}){
#        print "$cls\t$cls_s\t$anno\t$anno_s\t$anno_size{$anno}";
#        print "\n";
        my $pvalue = calculateStatistic(n11=>$anno_s, n1p=>$cls_s, np1=>$anno_size{$anno}, npp=>$M);
         # anno_term, anno_size, clsper, anno_total, backper, enrichment, pvalue
        push @tmp, [$anno, $anno_s, $anno_s/$cls_s, $anno_size{$anno}, $anno_size{$anno}/$M, $anno_s*$M/($cls_s*$anno_size{$anno}), $pvalue];
    #    push $resu{$cls}, [sort{$a[0] <=> $b[0]} @tmp];
    }
    @tmp = sort {$$a[6] <=> $$b[6]} @tmp;
    $resu{$cls} = [@tmp];
}

store \%resu, $store_file;
open(TMP, "> $anno_file") || die;
print TMP "ClsName\tClsSize\tAnno_term\tAnno_size\tClsPer\tAnno_total\tSeq_total\tBackPer\tEnrichment\tPvalue\n";
while(my ($cls, $info) = each %resu){
    foreach my $a (@{$info}){ #[$pvalue, $enrichment, $anno_s, $anno]
        print TMP join("\t",($cls, $cls_size{$cls}, $a->[0], $a->[1], $a->[2], $a->[3],
                        $M, $a->[4], $a->[5], $a->[6]))."\n";
#        print "$cls\t".join("\t",@{$a})."\n";
    }
#    print "$cls\t$#{$info}\n";
}
close(TMP)


