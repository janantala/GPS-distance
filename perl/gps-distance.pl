#!/usr/bin/perl

#:::  This routine calculates the distance between two points (given the
#:::  latitude/longitude of those points).
#:::  Definitions:                                                           :::
#:::    South latitudes are negative, east longitudes are positive           :::
#:::  Passed to function:                                                    :::
#:::    lat1, lon1 = Latitude and Longitude of point 1 (in decimal degrees)  :::
#:::    lat2, lon2 = Latitude and Longitude of point 2 (in decimal degrees)  :::
#:::    unit = the unit you desire for results                               :::
#:::           where: 'M' is statute miles (default)                         :::
#:::                  'K' is kilometers                                      :::
#:::                  'N' is nautical miles                                  :::

$pi = atan2(1,1) * 4;

sub distance {
	my ($lat1, $lon1, $lat2, $lon2, $unit) = @_;
	my $theta = $lon1 - $lon2;
	my $dist = sin(deg2rad($lat1)) * sin(deg2rad($lat2)) + cos(deg2rad($lat1)) * cos(deg2rad($lat2)) * cos(deg2rad($theta));
  $dist  = acos($dist);
  $dist = rad2deg($dist);
  $dist = $dist * 60 * 1.1515;
  if ($unit eq "K") {
  	$dist = $dist * 1.609344;
  } elsif ($unit eq "N") {
  	$dist = $dist * 0.8684;
		}
	return ($dist);
}

#:::  This function get the arccos function using arctan function   :::
sub acos {
	my ($rad) = @_;
	my $ret = atan2(sqrt(1 - $rad**2), $rad);
	return $ret;
}

#:::  This function converts decimal degrees to radians             :::
sub deg2rad {
	my ($deg) = @_;
	return ($deg * $pi / 180);
}

#:::  This function converts radians to decimal degrees             :::
sub rad2deg {
	my ($rad) = @_;
	return ($rad * 180 / $pi);
}

print distance(32.9697, -96.80322, 29.46786, -98.53506, "M") . " Miles\n";
print distance(32.9697, -96.80322, 29.46786, -98.53506, "K") . " Kilometers\n";
print distance(32.9697, -96.80322, 29.46786, -98.53506, "N") . " Nautical Miles\n";
