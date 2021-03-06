#!/usr/local/bin/perl

 use POSIX ();
 use POSIX qw(setsid);
 use POSIX qw(:errno_h :fcntl_h);

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval nanosleep
                      clock_gettime clock_getres clock_nanosleep clock
                      stat  );

my $start_time=[Time::HiRes::gettimeofday()];

#simulates one office with one window using mkSchedule, uSchedule, Radiance and EnergyPlus



# Office dimensions
$width=3; #X axis
$length=5; #Y axis
$height=2.7; #Z axis

#Window dimensions
$ledge=1; #height of the bottom part of the window (I did not know how to name it...?)
$wall_offset=0.1; #the window will be ($width-2*$wall_offset) wide
$ceiling_offset=0.1; #same deal, but with the ceiling. The window will be $height-$ledge-$ceiling_offset high



$sensor_spacing=0.5;
$sensor_height=0.8;

$epw_file="Colorado.epw"; #EnergyPlus weather file

$orientation=0; # 0 positions the window facing South, 90 facing East.

#Control.
$control="script.lua"; 
$n_positions=6; #number of positions of the shading system
$luminaire_power=30; #power of each luminaire (W)
$luminaire_file="lum.rad"; # 1 luminaria. Se asume que este archivo es .rad
$luminaire_spacing=1.5;
$lum_ceiling_offset=0.2; #espacio del ceiling al centro de la luminaria

$T_min=20; #thermostat setpoints.
$T_max=24; 

@control_sens=(0.5,0.5,0.8); #Fraccion del width, del length, altura (en m)

$cores=1; #for multicore processing (rcontrib, rtrace and rfluxmtx runs)
$bins=4; #number of bins of Reinhart's sky subdivition.

$N_people=0; #number of people in the office

$light_max=1000; # Range of acceptable illuminance (for plotting the illuminance results in Lightsolve Format)
$light_min=300; #


`mkdir -p scenes`;
`mkdir -p scenes/rad`;
`mkdir -p Workplanes`;
`mkdir -p DMX`;
`mkdir -p VMX`;
`mkdir -p LMX`;
`mkdir -p Schedules`;
`mkdir -p Results`;
`mkdir -p Results/andersen`;




$volume=$heigth*$width*$length;
$floor_area=$length*$width;
$window_area=($width-2*$wall_offset)*($height-$ceiling_offset-$ledge);
$wwr=$window_area/($height*$width);
print "\tWindow-to-wall Ratio: $wwr\n";


$lum_x=int($width/$luminaire_spacing); #this will create an array of luminaires. lum_x in the X axis and lum_y in the Y axis
$lum_y=int($length/$luminaire_spacing);
$n_luminaires=$lum_x*$lum_y;
$lighting_power=$n_luminaires*$luminaire_power;
$light_over_area=$lighting_power/$floor_area;
print "\tLuminaire power: $light_over_area W/m2\n\n";

$rest_x=($width-($lum_x-1)*$luminaire_spacing)/2;
$rest_y=($length-($lum_y-1)*$luminaire_spacing)/2;
$luminaire_height=$height-$lum_ceiling_offset;



my $t=Time::HiRes::tv_interval($start_time);
print "	Write Radiance model	$t\n";

$scene_file="Scene";
$z7=$height-$ceiling_offset;
$x8=$width-$wall_offset;
$z8=$z7;
$x9=$x8;


$rho_wall=0.7;
$rho_ceiling=0.7;
$rho_floor=0.7;

open (SCENE, ">>temp_rad");

print SCENE "void plastic wall\n0\n0\n5\t$rho_wall\t$rho_wall\t$rho_wall\t0\t0 \n\n";
print SCENE "void plastic floor\n0\n0\n5\t$rho_floor\t$rho_floor\t$rho_floor\t0\t0 \n\n";
print SCENE "void plastic ceiling\n0\n0\n5\t$rho_ceiling\t$rho_ceiling\t$rho_ceiling\t0\t0 \n\n";
	
print SCENE "wall polygon muro_trasero\n0\n0\n12\n\t0\t$length\t0\n\t$width\t$length\t0\n\t$width\t$length\t$height\n\t0\t$length\t$height\n\n";
print SCENE "wall polygon muro_derecho\n0\n0\n12\n\t$width\t0\t0\n\t$width\t$length\t0\n\t$width\t$length\t$height\n\t$width\t0\t$height\n\n";
print SCENE "wall polygon muro_izquierdo\n0\n0\n12\n\t0\t0\t0\n\t0\t$length\t0\n\t0\t$length\t$height\n\t0\t0\t$height\n\n";
print SCENE "ceiling polygon techo_geom\n0\n0\n12\n\t0\t0\t$height\n\t$width\t0\t$height\n\t$width\t$length\t$height\n\t0\t$length\t$height\n\n";
print SCENE "floor polygon piso_geom\n0\n0\n12\n\t0\t0\t0\n\t$width\t0\t0\n\t$width\t$length\t0\n\t0\t$length\t0\n\n";

print SCENE "wall polygon muro_ventana
0
0
30
	0	0	0
	$width	0	0
	$width	0	$height
	0	0	$height
	0	0	0
	$wall_offset	0	$ledge
	$wall_offset	0	$z7
	$x8	0	$z8
	$x9	0	$ledge
	$wall_offset	0	$ledge

";


close (SCENE); 	

`xform -t $rest_x $rest_y $luminaire_height -a $lum_x -t $luminaire_spacing 0 0 -a $lum_y -t 0 $luminaire_spacing 0 $luminaire_file >> temp_rad`;

open (SCENEFILE, "> scenes/rad/$scene_file.rad");


close (SCENEFILE);
`xform -rz $orientation temp_rad >> scenes/rad/$scene_file.rad`;
`rm -f temp_rad`;


$win_file="window";
open (WIN_SCENE, "> win_rad");
print WIN_SCENE "#\@rfluxmtx h=kf u=Z \nvoid glass black\n0\n0\n3\t0\t0\t0\n\n \nblack\tpolygon\twindow\n0\n0\n12\t$wall_offset\t0\t$ledge\n\t$wall_offset\t0\t$z7\n\t$x8\t0\t$z7\n\t$x8\t0\t$ledge\n";
close (WIN_SCENE);
`xform -rz $orientation win_rad > scenes/rad/$win_file.rad`;
`rm -f win_rad`;
my $t=Time::HiRes::tv_interval($start_time);
print "	Write window file	$t\n";




$control_sens="Workplanes/control.pts";


open (SENS, ">".$control_sens);
@tmp_array=@control_sens;

while ($#tmp_array >= 0) {
	$x=@tmp_array[0];
	$y=@tmp_array[1];
	$z=@tmp_array[2];
	
	$a=$x*$width*cos($orientation*3.141592654/180)-$y*$length*sin($orientation*3.141592654/180);
	$b=$x*$width*sin($orientation*3.141592654/180)+$y*$length*cos($orientation*3.141592654/180);		

	print SENS "$a $b $z 0 0 1\n";
	shift @tmp_array;
	shift @tmp_array;
	shift @tmp_array;
}

close(SENS);
my $t=Time::HiRes::tv_interval($start_time);
print "	Write control sensors	$t\n";

$sens_x=int($width/$sensor_spacing);
$sens_y=int($length/$sensor_spacing);

$rest_x=($width-($sens_x-1)*$sensor_spacing)/2;
$rest_y=($length-($sens_y-1)*$sensor_spacing)/2;


$wp="Workplanes/WP_$length-$width-$orientation.pts";
if (not -e $wp) {
	open (SENS, ">>".$wp);
	
	for ($i=$rest_x; $i<$width; $i+=$sensor_spacing){
		for ($j=$rest_y; $j<$length; $j+=$sensor_spacing){			
			$a=$i*cos($orientation*3.141592654/180)-$j*sin($orientation*3.141592654/180);
			$b=$i*sin($orientation*3.141592654/180)+$j*cos($orientation*3.141592654/180);		
			print SENS "$a $b $sensor_height 0 0 1\n";
		}
	}
	
	close SENS;
}
my $t=Time::HiRes::tv_interval($start_time);
print "	Write Workplane sensors	$t\n";



# Se asume que no hay otras fuentes de luz.
$lum_wp="LMX/LUM_WP.lmx";

`oconv scenes/rad/$scene_file.rad scenes/rad/$win_file.rad > octree`;
`cat $wp | rtrace -n $cores -I -af af -ab 6 -ad 2048 -aa 0.2 octree > $lum_wp`;
`rm octree af`;
my $t=Time::HiRes::tv_interval($start_time);	
print "	Luminaire contribution to workplane	$t\n";




$lum_c="LMX/LUM_Control.lmx";

`oconv scenes/rad/$scene_file.rad scenes/rad/$win_file.rad > octree`;
`cat $control_sens | rtrace -n $cores -I -af af -ab 6 -ad 2048 -aa 0.2 octree > $lum_c`;
`rm octree af`;
my $t=Time::HiRes::tv_interval($start_time);
print "	Luminaire contribution to control sensors	$t\n";


## Scene file
$dmx_file="DMX/D_$scene_file.dmx";
$or_x=sin($orientation*3.141592654/180);
$or_y=-cos($orientation*3.141592654/180);
$or_z=0;

`echo "#\@rfluxmtx h=u u=Y\nvoid glow ground_glow\n0\n0\n4 1 1 1 0\n\nground_glow source ground\n0\n0\n4 0 0 -1 180\n\n#\@rfluxmtx h=r$bins u=Y\n\nvoid glow skymat\n0\n0\n4 1 1 1 0\n\nskymat source sky\n0\n0\n4 0 0 1  180\n\n" > white_sky`;
#`oconv  scenes/rad/$scene_file.rad white_sky > octree`;
#`genklemsamp -c 1000 -vd $or_x $or_y $or_z scenes/rad/$win_file.rad | rcontrib -n $cores -c 1000 -e MF:$bins -f reinhart.cal -b rbin -bn Nrbins -m skymat octree > $dmx_file`;
`rm -f white_sky octree`;
`rfluxmtx scenes/rad/$win_file.rad white_sky > $dmx_file`;
my $t=Time::HiRes::tv_interval($start_time);
print "	Calculate Daylight matrix	$t\n";





$vmx_file_WP="VMX/WP_$scene_file.vmx";
$or_x=-$or_x;
$or_y=-$or_y;

`echo "void glow winmat\n0\n0\n4 1 1 1 0\n\n void glass black_glass\n0\n0\n3 0 0 0\n\n" > winfile`;
`xform -m winmat scenes/rad/$win_file.rad >> winfile`;
#`oconv winmat winfile scenes/rad/$scene_file.rad > octree`;
#`cat $wp | rcontrib -f klems_int.cal -b 'kbin($or_x,$or_y,0,0,0,1)' -bn Nkbins -m winmat -I+ -n $cores -ab 6 -ad 2048 -aa 0.2 octree > $vmx_file_WP`;
`rfluxmtx -I+ -n $cores -ab 9 -ad 16384 -lw 6.1e-5 <$wp - winfile scenes/rad/$scene_file.rad > $vmx_file_WP`;
`rm -f winmat winfile octree af`;

my $t=Time::HiRes::tv_interval($start_time);
print "	Calculate View matrix over workplane	$t\n";


##
$vmx_file_C="VMX/C_$scene_file.vmx";
`echo "void glow winmat\n0\n0\n4 1 1 1 0\n\n void glass black_glass\n0\n0\n3 0 0 0\n\n" > winfile`;
`xform -m winmat scenes/rad/$win_file.rad >> winfile`;
#`oconv winmat winfile scenes/rad/$scene_file.rad > octree`;
#`cat $control_sens | rcontrib -f klems_int.cal -b 'kbin($or_x,$or_y,0,0,0,1)' -bn Nkbins -m winmat -I+ -n $cores -ab 6 -ad 2048 -aa 0.2 octree > $vmx_file_C`;
`rfluxmtx -I+ -n $cores -ab 9 -ad 16384 -lw 6.1e-5 <$control_sens - winfile scenes/rad/$scene_file.rad > $vmx_file_C`;
`rm -f winmat winfile octree af`;

my $t=Time::HiRes::tv_interval($start_time);
print "	Calculate View matrix over control sensors	$t\n";



##
$schedule="Schedules/schedule.txt";
`cat $lum_c > LMX-control-1.lmx`; 
`cat $vmx_file_C > WindowSet_1-control.vmx`;
`cat $dmx_file > WindowSet_1.dmx`;
`./mkSchedule -o AAA -f $epw_file -m $bins -w 1 -x 1 -l 1 -L LMX-control-%d.lmx -V WindowSet_%d-control.vmx -D WindowSet_%d.dmx -T 6 BSDF/blind%d-klems.xml -u $control > $schedule`;
`rm LMX-control-1.lmx WindowSet_1-control.vmx WindowSet_1.dmx`;

my $t=Time::HiRes::tv_interval($start_time);
print "	run mkSchedule	$t\n";


$light_file="Results/results.txt";

$x=$sens_y*$sens_x;

`cat $lum_wp > LMX-workplane-1.lmx`; 
`cat $vmx_file_WP > WindowSet_1-WP.vmx`;
`cat $dmx_file > WindowSet_1.dmx`;
`./uSchedule -n 8760 -u $schedule -V WindowSet_%d-WP.vmx -L LMX-workplane-%d.lmx -x $x -m $bins > $light_file`;
`rm LMX-workplane-1.lmx WindowSet_1-WP.vmx WindowSet_1.dmx`;

my $t=Time::HiRes::tv_interval($start_time);
print "	Lighting simulation (uSchedule)	$t\n";


## Data analysis: create a Temporal Map in Lightsolve FOrmat (red is "too much", green is "ok" and blue is "too low").

$andersen_file="Results/andersen/$ciudad-$control-$light_max-$light_min.hdr";


	### Mapa Temporal Andersen
	$hour=8;
	
	$x=24*$hour;
	$y=365;

	open(FILE,$light_file); 
	my @lines=<FILE>;
	close(FILE);


	$header="#?RADIANCE\n\n-Y $y +X $x\n";

	open (MYFILE, '>', $andersen_file);
	print MYFILE $header;

	foreach $ln (@lines){
		
		my @andersen=split(",",$ln);
		$under=0;
		$good=0;
		$over=0;
		$n=0;
		
		foreach $set (@andersen){
			$n++;
			if($set > $light_max){
				$over++;
			}elsif($set < $light_min){
				$under++;
			}else{
				$good++;
			}
		}
		
		$good=sprintf("%f",$good/$n);
		$under=sprintf("%f",$under/$n);
		$over=sprintf("%f",$over/$n);

		my ($rr, $gr, $br, $e)=rgbe($over,$good,$under);

		for my $i (1..$hour){
			print MYFILE pack('C',$rr);
			print MYFILE pack('C',$gr);
			print MYFILE pack('C',$br);
			print MYFILE pack('C',$e);
		}
		
	}

	close (MYFILE);
	
	`ra_tiff $andersen_file andersen.tiff`;





								








sub rgbe{
	my ($r, $g, $b)=@_;

	$rr;
	$gr;
	$br;
	$e;

	$d=-10000;

	if($r>$d){$d=$r;}
	if($g>$d){$d=$g;}
	if($b>$d){$d=$b;}

	if($d <= 1e-32){
		$rr=0;
		$gr=0;
		$br=0;
		$e=0;
	}else{
		($m, $e) = POSIX::frexp($d);
		$d=$m*255.999/$d;

		if($r>0){
			$rr=$r*$d;
			$rr=int($rr + $rr/abs($rr*2));
		}else{
			$rr=0;
		}

		if($g>0){
					$gr=$g*$d;
					$gr=int($gr + $gr/abs($gr*2));
			}else{
					$gr=0;
			}
	
		if($b>0){
					$br=$b*$d;
					$br=int($br + $br/abs($br*2));
			}else{
					$br=0;
			}

		$e=$e+128;
		$e=int($e + $e/abs($e*2));
	}
	

	return ($rr, $gr, $br, $e);
	
}
