#!usr/bin/perl -w
#
#liulin@frasergen.com
#
open IN, "$ARGV[0]";
my $start=1;
my $end=$ARGV[1];
my $chr=$ARGV[2];
while (my $line = <IN>){
	chomp $line;
	next if ($line =~ /^#/);
	next if ($line =~ /^header/);
	my @a = split(/\s+/,$line);
	$a[1] = $a[1] -1;
	print "$chr\t$start\t$a[1]\n";
	$start=$a[2];
}
close IN;
print "$chr\t$start\t$end\n";
