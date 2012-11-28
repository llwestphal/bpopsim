###
# Pod Documentation
###

=head1 NAME

KS_test_2d

=head1 SYNOPSIS

Performs the 2-dimensional Kolmogorov-Smirnov test as described in [Fasano, G.
and Franceschini, A. A multidimensional version of the Kolmogorov-Smirnov test. 
Monthly Notes of the Royal Astronomical Society (1987). 255: 155-170.]

=head1 AUTHOR

Jeffrey Barrick <jbarrick@msu.edu>

=head1 COPYRIGHT

Copyright (C) 2008-2009.

=cut

###
# End Pod Documentation
###

package KS_test_2d;
use strict;

require Exporter;
our @ISA = qw( Exporter );
our @EXPORT = qw( bf_KS_test_2d z_value_at_significance_level d_value_2d_to_probability);

use Data::Dumper;

=head2 KS_test_2d

 Title   : bf_KS_test_2d
 Usage   : $sig = KS_test_2d( \@list_1, \@list_2 );
 Function: 
 Returns : brute force two dimensional Kolmogorov-Smirnov significance level p-value
 
=cut

#brute force version.
sub bf_KS_test_2d
{
	my ($list_1_ref, $list_2_ref) = @_;
	
	#tell us which each is from
	foreach my $item (@$list_1_ref)
	{
		$item->{type} = 1;
	}
	foreach my $item (@$list_2_ref)
	{
		$item->{type} = 2;
	}
	my @list = (@$list_1_ref, @$list_2_ref);
	
	##calculate the correlation coefficient for the combined list
	my $cc = pearson_correlation_coefficient(\@list);
	my $cc1 = pearson_correlation_coefficient($list_1_ref);
	my $cc2 = pearson_correlation_coefficient($list_2_ref);
	
	print "Correlation coefficient: list_1 = $cc1, list_2 = $cc2, combined = $cc\n";
	
	#now we cycle through each point, choosing it as the location to draw the quadrants
	my $d;
	foreach(my $i=0; $i<scalar @list; $i++)
	{
		my $center = $list[$i];
		
		#count up statistics for each quadrant
		my @d_vals;
		for (my $k=0; $k<4; $k++)
		{
			$d_vals[$k]->{1} = 0;
			$d_vals[$k]->{2} = 0;
		}		 
		 		 
		foreach(my $j=0; $j<scalar @list; $j++)
		{
			next if ($i == $j);
			my $test = $list[$j];
			
			my $q;
			if ($test->{y} > $center->{y})
			{
				$q = ($test->{x} > $center->{x}) ? 0 : 1;
			}
			else
			{
				$q = ($test->{x} > $center->{x}) ? 2 : 3;
			}
			
			$d_vals[$q]->{$test->{type}}++;
		}
		
		#add center point to each quadrant
		for (my $k=0; $k<4; $k++)
		{
			$d_vals[4+$k]->{1} = $d_vals[$k]->{1};
			$d_vals[4+$k]->{2} = $d_vals[$k]->{2};
			$d_vals[4+$k]->{$center->{type}}++;
		}
		
		foreach my $d_val (@d_vals)
		{
			$d_val = abs( ($d_val->{1} / scalar @$list_1_ref) - ($d_val->{2} / scalar @$list_2_ref));
			$d = $d_val if (!defined $d || ($d_val) > $d);
		}
	}
	
	print " D = $d\n";
	
	return d_value_2d_to_probability($d, scalar @$list_1_ref, scalar @$list_2_ref, $cc);
}

sub max_from_list
{
	my $max = pop @_;
	foreach my $i (@_)
	{
		$max = $i if ($i > $max);
	}
	return $max;
}

sub pearson_correlation_coefficient
{
	my ($list_ref) = @_;
	my $n = scalar @$list_ref;

	## calculate means
	my $mean_x = 0;
	my $mean_y = 0;
	foreach my $i (@$list_ref)
	{
		$mean_x += $i->{x};
		$mean_y += $i->{y};
	}
	$mean_x /= $n;
	$mean_y /= $n;
	
	## calculate sample standard deviatinos
	my $sd_x = 0;
	my $sd_y = 0;
	foreach my $i (@$list_ref)
	{
		$sd_x += ($i->{x} - $mean_x)**2;
		$sd_y += ($i->{y} - $mean_y)**2;
	}
	$sd_x = sqrt($sd_x/($n-1));
	$sd_y = sqrt($sd_y/($n-1));
	
	## calculate correlation coefficient
	my $cc = 0;
	foreach my $i (@$list_ref)
	{
		$cc += (($i->{x} - $mean_x) / $sd_x) * (($i->{y} - $mean_y) / $sd_y);
	}
	$cc /= $n-1;
	
	return $cc;
}

#binary tree version (not completed).
sub bt_KS_test_2d
{
	my ($list_1_ref, $list_2_ref) = @_;
	my @list;
	my $sl;
	my $z;

	### merge lists and sort
	@list = sort { $a->{y} <=>  $b->{y} } (@$list_1_ref, @$list_2_ref);
	my $root;
	while (@list)
	{
		my $new_node = pop @list;
		
		#add new node to a last child 
		$root = addnode($root, $new_node);
		
		#propagate scores up
		propagate_upward($new_node);
	}
	
	
	
	### now sort the other way
	@list = sort { -($a->{y} <=>  $b->{y}) } (@$list_1_ref, @$list_2_ref);
	
	return ($sl, $z);
}

sub add_node
{
	my ($root, $node) = @_;
	
	return $node if (!defined $root) ;

	my $search = $node;
	while ($search)
	{
		if ($search->{x} < $node->{x})
		{
			if (!defined $search->{left})
			{
				$search->{left} = $node;
				$node->{parent} = $search;
				return $root;
			}
			else
			{
				$search = $search->{left};
			}
		}
		else # if ($search->{x} >= $node->{x})
		{
			if (!defined $search->{right})
			{
				$search->{right} = $node;
				$node->{parent} = $search;
				return $root;
			}
			else
			{
				$search = $search->{right};
			}
		}
	}
	
	die "Could not add node!" if (!defined $search);
}

sub propagate_upward
{
	my ($start_node) = @_;
	
	for (my $on_node = $start_node; $on_node; $on_node = $on_node->{parent_node})
	{
	}
}


sub z_value_1d_to_probability
{
	my ($z) = @_;
	
	my $pr_gt_z = 0;
	for (my $i=1; $i<1000; $i++)
	{
		$pr_gt_z += (-1)**($i-1) * exp(-2 * $i**2 * $z);
		#print "$pr_gt_z\n";
	}
	$pr_gt_z *= 2;
	return $pr_gt_z;
}

sub d_value_2d_to_probability
{
	my ($d, $n1, $n2, $cc) = @_;
	
	my $n = $n1*$n2 / ($n1 + $n2);
	my $val = $d * sqrt($n) / (1 + sqrt(1-$cc**2) * (0.25 - 0.75/sqrt($n)));
	
	print STDERR "$d, $n1, $n2, $cc // $n, $val\n";
	return z_value_1d_to_probability($val)
}

####
#### BEGIN UNUSED CODE
####

#internal function, creates a table of Z values
sub create_z_table
{
	my ($n, $cc, $precision) = @_;
	$precision = 0.1 if (!defined $precision);	

	my $table;
	$table->{precision} = $precision;
	for (my $sl = 0; $sl < 100; $sl += $precision)
	{
		push @{$table->{list}}, z_value_at_significance_level($n, $cc, $sl);	
	}
	
	return $table;
}

#internal function, looks up given Z value in table
sub z_value_to_significance_level
{
	my ($z, $table) = @_;
}

# Fasano, G. and Franceschini, A. A multidimensional version of the Kolmogorov-
# Smirnov test. Monthly Notes of the Royal Astronomical Society (1987). 255: 155-170.
our $coefficients = [

# i=0, j=0
	[
		[ 0.7107, 0.009853, 0.008561, 0.002800 ],
# i=0, j=1
		[ -0.5126, -0.2539, -0.03745 ],
# i=0, j=2
		[ -0.3298, -0.01408 ],
# i=0, j=3
		[ -0.07693 ]
	],
# i=1, j=0
	[
		[ -0.06124, -0.01918, 0.01064 ],
# i=1, j=1
		[ -0.05042, 0.03670 ],
# i=1, j=2
		[ -0.0007750 ],
	],
# i=2, j=0
	[
		[ 0.7775, -0.2322 ],
# i=2, j=1
		[ -0.1273 ]
	],
# i=3, j=0
	[
		[ -0.9915 ] 
	]
];

#internal function, calculates the Z value corresponding to the significance level
## NOTE: THIS DOES NOT WORK CORRECTLY!?? 
sub ff_z_value_at_significance_level
{
	my ($n, $cc, $sl) = @_;

	my $s = log($n)/log(10) + 1.074;
	my $u = -log((1.268 - $cc)/1.41)/log(10);
	my $v = log(100-$sl)/log(10) - 2;	
	my $w = $s * (4.989 / (4.989 - $s))**-0.09418;

#	my $s = log($n)/log(10) + 1.074;
#	my $u = -log((1.268 - $cc)/1.41)/log(10);
#	my $v = log(100-$sl)/log(10) - 2;
#	my $w = $s * exp (-0.09418 *  log(4.989 / (4.989 - $s)));

	print "n = $n, cc = $cc, sl = $sl\n";
	print "u = $u, v = $v, w = $w, s = $s\n";

	#print Dumper($coefficients);

	my $z = 0;
	for (my $i=0; $i<=3; $i++)
	{
		for (my $j=0; $j<=3-$i; $j++)
		{
			for (my $k=0; $k<=3-$i-$j; $k++)
			{
				print "$i $j $k $coefficients->[$i]->[$j]->[$k]\n";				
				print "   u^i = " . $u**$i . "\n";		
				print "   v^j = " . $v**$j . "\n";		
				print "   w^k = " . $w**$k . "\n";		

				my $inc = $coefficients->[$i]->[$j]->[$k] * $u**$i * $v**$j * $w**$k;
				print "  INC = $inc\n";
				$z += $inc;
				print "  Z = $z\n";
			}
		}
	}
	return $z;
}

return 1;
