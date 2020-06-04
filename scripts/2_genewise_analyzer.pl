#!/usr/bin/perl
use strict;
use warnings;
my @infiles=glob("*_sitelk_table.txt");
foreach my $infile (@infiles) {
  next unless $infile=~m/(?<branch>\w+)_Gene_(?<id>\w+)_sitelk_table.txt/gi;
  my $gene_name=$+{id};
  my $branch_name=$+{branch};
  print "$branch_name\t$gene_name\n";

  my $t1_lk=0;
  my $t2_lk=0;
  my @tree1s;
  my @tree2s;
  open (MATCHIN, $infile) or die "Can't open the file $infile: $!";
  while (my $content=<MATCHIN>) { 
	chomp $content;
    my @elements=split(/\t+/,$content);
    $t1_lk=$t1_lk+$elements[5];
	$t2_lk=$t2_lk+$elements[6];
	push @tree1s,$elements[3] if $elements[3]eq "tr1";
	push @tree2s,$elements[3] if $elements[3]eq "tr2";
  }
  close MATCHIN;
  my $tree1s_no=@tree1s;
  my $tree2s_no=@tree2s;
  my $out= "${branch_name}_sitwise_statistics.txt";
  open (OUT,">>", $out) or die "Can't write the file $out: $!";
  print OUT "$gene_name\ttr1\t$tree1s_no\ttr2\t$tree2s_no\n";
  close OUT;
  my $outputs= "${branch_name}_genewise_lnL.txt";
  open (OUTPUT,">>", $outputs) or die "Can't write the file $outputs: $!";
  if ($t1_lk>$t2_lk) {print OUTPUT "$gene_name\tt1\tt1\t$t1_lk\tt2\t$t2_lk\n";}
  elsif($t2_lk>$t1_lk) {print OUTPUT "$gene_name\tt2\tt1\t$t1_lk\tt2\t$t2_lk\n";}
  close OUTPUT;
 }
