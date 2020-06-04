#!/usr/bin/perl
use strict;
use warnings;
 my @infiles=glob("RAxML_perSiteLLs.*");
 foreach my $infile (@infiles) {
  next unless $infile=~m/RAxML_perSiteLLs.(?<id>[a-z]+)_(?<branch>\w+)_site_lk/gi;
  my $id_name=$+{id};
  my $branch_name=$+{branch};

  my %gene;
  my %gene_position;
  my $filename="${id_name}_combining_orders.txt";
  open (MATCHIN, $filename) or die "Can't open the file $filename: $!";
  while (my $content=<MATCHIN>) { 
	chomp $content;
 	$content=~s/.*,\s+(.*)/$1/gi;
	$content=~s/-/ /;
	print "$content\n";
    my @elements=split(/\s+/,$content);
    my $gene_name=$elements[0];
    my $gene_start=$elements[2];
    my $gene_end=$elements[3];
	 $gene_position{$gene_name}=$gene_start-1;
    my @gene_sites=($gene_start..$gene_end);
    foreach  my $gene_site( @gene_sites) { $gene{$gene_site}=$gene_name;}		 
  }
 close MATCHIN;





  my %site_hash;
  my %tree_hash;
  my %sites_hash;
  my %trees_hash;
  my @trees;
  my @all_lns;
  open (IN, $infile) or die "Can't open the file $infile: $!";
  while (my $content=<IN>)  {
	 chomp $content;
	 $content=~s/\t+/ /gi;
	 next unless $content=~m/^(?<tre>t\w+\d+)\s+(?<lns>.*)/gi;
	 my $tree_name=$+{tre};
	 push  @trees, $tree_name;
	 my $lns_content=$+{lns};
     my @sites_lns=split(/\s+/,$lns_content); 
	 my $site=0;
     foreach my $site_ln (@sites_lns) {
	   print "reading your data...\n";
	   $site++;
	   my $i=$site-1;
	   push @all_lns, $site_ln;
	   ## recognize putative site with lighest value ##
	   if (exists $site_hash{$site}) {
        if ($site_hash{$site}<$site_ln) {$site_hash{$site}=$sites_lns[$i];$tree_hash{$site}=$tree_name;}
	      }else{ $site_hash{$site}=$sites_lns[$i];$tree_hash{$site}=$tree_name;}

       ## joint all vaules for each site ##
	   if (exists $sites_hash{$site}) {
		  $sites_hash{$site}=join("\t",$sites_hash{$site},$sites_lns[$i]);$trees_hash{$site}=join(" ",$trees_hash{$site},$tree_name);
	    }else{ $sites_hash{$site}=$sites_lns[$i];$trees_hash{$site}=$tree_name;}
	 }  
  }
 close IN;




 
   foreach my $decisive_site (sort {$a <=>$b} keys %tree_hash) {
	   my $gene_name=$gene{$decisive_site};
	   my $site_rela_pos=$decisive_site-$gene_position{$gene_name};
	   my $output2= "${branch_name}_Gene_${gene_name}_sitelk_table.txt";
      open (OUTPUT2,">>", $output2) or die "Can't open the file $output2: $!";
	  print OUTPUT2 "$gene_name\t$decisive_site\t$site_rela_pos\t$tree_hash{$decisive_site}\t$trees_hash{$decisive_site}\t$sites_hash{$decisive_site}\n"; 
	  close OUTPUT2;
   }
   


  
  print "done\n";
}


