#!perl -w

use 5.014;
use strict;
use autodie;

my $se_1_ref = "./SRR065390_1_sub.fastq";
my $se_2_ref = "./SRR065390_2_sub.fastq";
my $se_1_act = "./decom.SRR065390_1_sub.fastq";
my $se_2_act = "./decom.SRR065390_2_sub.fastq";

#my $se_1_act = "./rd_1.jnk.align";
#my $se_2_act = "./rd_2.jnk.align";
#my $se_1_act = "./rd_1.jnk.compact";
#my $se_2_act = "./rd_2.jnk.compact";

my $summary = "./summary.txt";
my $sort_ends = 1;
my $sort_reads = 1;


my $one_read;
my @ref_pe_reads;
my @act_pe_reads;
my %ref_se_reads;
my @act_se_reads;


my $se_1_ref_rl = -1;
my $se_1_ref_rd_cnt = 0;
open SE_1_REF, '<', $se_1_ref;
while (my $a_line = <SE_1_REF>) {
    if ($a_line =~ m/\@/) {
        #FASTQ format
        chomp($one_read = <SE_1_REF>); #Read line
        $a_line = <SE_1_REF>; #Comment line
        die("Incorrect file format") unless ($a_line =~ m/^\+/);
        $a_line = <SE_1_REF>; #Quality string
    } elsif ($a_line =~ m/\>/) {
        #FASTA format
        chomp($one_read = <SE_1_REF>); #Read line
    } else {
        die("Should not be here");
    }
    
    if ($se_1_ref_rl == -1) {
        $se_1_ref_rl = length($one_read);
        say "1st read in $se_1_ref - $one_read";
    }
    die("Unequal read lengths") unless ($se_1_ref_rl == length($one_read));
    $ref_pe_reads[$se_1_ref_rd_cnt] = "\U$one_read"; #Read from se_1_ref
    $ref_se_reads{"\U$one_read"} = $se_1_ref_rd_cnt; #Read from se_1_ref
    $se_1_ref_rd_cnt++;
}
close SE_1_REF;
say "Processed $se_1_ref_rd_cnt reads from $se_1_ref";


my $se_2_ref_rl = -1;
my $se_2_ref_rd_cnt = 0;
open SE_2_REF, '<', $se_2_ref;
while (my $a_line = <SE_2_REF>) {
    if ($a_line =~ m/\@/) {
        #FASTQ format
        chomp($one_read = <SE_2_REF>); #Read line
        $a_line = <SE_2_REF>; #Comment line
        die("Incorrect file format") unless ($a_line =~ m/^\+/);
        $a_line = <SE_2_REF>; #Quality string
    } elsif ($a_line =~ m/\>/) {
        #FASTA format
        chomp($one_read = <SE_2_REF>); #Read line
    } else {
        die("Should not be here");
    }
    
    if ($se_2_ref_rl == -1) {
        $se_2_ref_rl = length($one_read);
        say "1st read in $se_2_ref - $one_read";
    }
    die("Unequal read lengths") unless ($se_2_ref_rl == length($one_read));
    $ref_pe_reads[$se_2_ref_rd_cnt] = $ref_pe_reads[$se_2_ref_rd_cnt] . " " . "\U$one_read"; #Read from se_2_ref
    $ref_se_reads{"\U$one_read"} = $se_2_ref_rd_cnt + $se_1_ref_rd_cnt; #Read from se_1_ref
    $se_2_ref_rd_cnt++;
}
close SE_2_REF;
say "Processed $se_2_ref_rd_cnt reads from $se_2_ref";

die("Unequal read lengths") unless ($se_1_ref_rl == $se_2_ref_rl);
#say $ref_pe_reads[ 0];
##say $ref_pe_reads[ 1];
##say $ref_pe_reads[-2];
#say $ref_pe_reads[-1];


my $se_1_act_rl = -1;
my $se_1_act_rd_cnt = 0;
open SE_1_ACT, '<', $se_1_act;
while (my $a_line = <SE_1_ACT>) {
    if ($a_line =~ m/\@/) {
        #FASTQ format
        chomp($one_read = <SE_1_ACT>); #Read line
        $a_line = <SE_1_ACT>; #Comment line
        die("Incorrect file format") unless ($a_line =~ m/^\+/);
        $a_line = <SE_1_ACT>; #Quality string
    } elsif ($a_line =~ m/\>/) {
        #FASTA format
        chomp($one_read = <SE_1_ACT>); #Read line
    } else {
        die("Should not be here");
    }
    
    if ($se_1_act_rl == -1) {
        $se_1_act_rl = length($one_read);
        say "1st read in $se_1_act - $one_read";
    }
    die("Unequal read lengths") unless ($se_1_act_rl == length($one_read));
    $act_pe_reads[$se_1_act_rd_cnt] = "\U$one_read"; #Read from se_1_act
    $act_se_reads[$se_1_act_rd_cnt] = "\U$one_read"; #Read from se_1_act
    $se_1_act_rd_cnt++;
}
close SE_1_ACT;
say "Processed $se_1_act_rd_cnt reads from $se_1_act";


my $se_2_act_rl = -1;
my $se_2_act_rd_cnt = 0;
open SE_2_ACT, '<', $se_2_act;
while (my $a_line = <SE_2_ACT>) {
    if ($a_line =~ m/\@/) {
        #FASTQ format
        chomp($one_read = <SE_2_ACT>); #Read line
        $a_line = <SE_2_ACT>; #Comment line
        die("Incorrect file format") unless ($a_line =~ m/^\+/);
        $a_line = <SE_2_ACT>; #Quality string
    } elsif ($a_line =~ m/\>/) {
        #FASTA format
        chomp($one_read = <SE_2_ACT>); #Read line
    } else {
        die("Should not be here");
    }
    
    if ($se_2_act_rl == -1) {
        $se_2_act_rl = length($one_read);
        say "1st read in $se_2_act - $one_read";
    }
    die("Unequal read lengths") unless ($se_2_act_rl == length($one_read));
    $act_pe_reads[$se_2_act_rd_cnt] = $act_pe_reads[$se_2_act_rd_cnt] . " " . "\U$one_read"; #Read from se_2_act
    $act_se_reads[$se_2_act_rd_cnt + $se_1_act_rd_cnt] = "\U$one_read"; #Read from se_1_act
    $se_2_act_rd_cnt++;
}
close SE_2_ACT;
say "Processed $se_2_act_rd_cnt reads from $se_2_act";

die("Unequal read lengths") unless ($se_1_act_rl == $se_2_act_rl);
#say $act_pe_reads[ 0];
##say $act_pe_reads[ 1];
##say $act_pe_reads[-2];
#say $act_pe_reads[-1];
die("Unequal read lengths") unless ($se_1_ref_rl == $se_2_act_rl);


sub by_special {
    my @var1 = split(/ /, $a);
    my @var2 = split(/ /, $b);
    if ($var1[0] gt $var1[1]) {
        ($var1[0], $var1[1]) = ($var1[1], $var1[0]);
    }
    if ($var2[0] gt $var2[1]) {
        ($var2[0], $var2[1]) = ($var2[1], $var2[0]);
    }
    
    if ($var1[0] lt $var2[0]) {
        return -1;
    } elsif ($var1[0] gt $var2[0]) {
        return 1;
    } else {
        if ($var1[1] lt $var2[1]) {
            return -1;
        } elsif ($var1[1] gt $var2[1]) {
            return 1;
        } else {
            $var1[2] <=> $var2[2];
        }
    }
}

if ($sort_reads) {
    my $var2;
    $var2 = 1;
    foreach my $var1 (@ref_pe_reads) {
        $var1 = $var1 . " " . $var2;
        $var2++;
    }
    $var2 = 1;
    foreach my $var1 (@act_pe_reads) {
        $var1 = $var1 . " " . $var2;
        $var2++;
    }
#    say $ref_pe_reads[-2];
#    say $act_pe_reads[-1];
    
    @ref_pe_reads = sort by_special @ref_pe_reads;
    @act_pe_reads = sort by_special @act_pe_reads;
#    say $ref_pe_reads[-2];
#    say $act_pe_reads[-1];
    
##    @ref_se_reads = sort @ref_se_reads;
##    @act_se_reads = sort @act_se_reads;
}


my $mat_pe_reads = 0;
my $mis_pe_reads = 0;
my $mat_se_reads = 0;
my $mis_se_reads = 0;
open SUMMARY, '>', $summary;
foreach my $var0 (0 .. $#ref_pe_reads) {
    my @var1 = split(/ /, $ref_pe_reads[$var0]);
    my @var2 = split(/ /, $act_pe_reads[$var0]);
    if ($sort_ends) {
        if ($var1[0] gt $var1[1]) {
            ($var1[0], $var1[1]) = ($var1[1], $var1[0]);
        }
        if ($var2[0] gt $var2[1]) {
            ($var2[0], $var2[1]) = ($var2[1], $var2[0]);
        }
    }
    
    if (($var1[0] eq $var2[0]) and ($var1[1] eq $var2[1])) {
        $mat_pe_reads++;
    } else {
#        if ($var1[0] eq $var2[0]) {
#            $mat_se_reads += 1;
#            $mis_se_reads += 1;
#        } elsif ($var1[1] eq $var2[1]) {
#            $mat_se_reads += 1;
#            $mis_se_reads += 1;
#        } else {
#            $mis_se_reads += 2;
#        }
        $mis_pe_reads++;
        say SUMMARY join(" ", @var1);
        say SUMMARY join(" ", @var2), "\n";
    }
}
foreach my $var0 (@act_se_reads) {
    if (exists $ref_se_reads{$var0}) {
        $mat_se_reads += 1;
    } else {
        $mis_se_reads += 1;
    }
}
close SUMMARY;

say "#" x 48;
say "Total PE read count: ", $mat_pe_reads + $mis_pe_reads;
say "Match PE read count: ", $mat_pe_reads;
say "Diff. PE read count: ", $mis_pe_reads;
say "Total SE read count: ", $mat_se_reads + $mis_se_reads;
say "Match SE read count: ", $mat_se_reads;
say "Diff. SE read count: ", $mis_se_reads;
say "#" x 48;



































