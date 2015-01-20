-- Germán Molina germolinal@gmail.com

-- This is an example of a Lua script that can be used as a control algorithm.


--The functions you can use in the Lua scripts to retrieve information are:

--hour() --> returns the actual time (hour)
--irradiance(x,y,z) --> returns the solar irradiance over a surface with normal [x,y,z]. It is according with the book Solar 
--Engineering of THermal Processes... not the same as EnergyPlus, but close. Assumes there are no obstructions.
--night(n) --> receives the nighttime illuminance over a sensor (the nth sensor). THat is, the illuminance of the sensor with all the lights on, and no sun.
--update() --> updates the illuminance values, useful for dimming the luminaires after changing the 
--sensor(n) --> returns the illuminance value on the nth sensor at the moment.
--ext_temp() --> returns the exterior temperature
--lum(n) -> returns the power fraction of the nth luminaire set (0 is off, 1 is fully on)
--win(n) --> returns the position (as an integer) of the nth window set.
--is_day() --> returns 1 if it is day, 0 if it is night.
--azimuth() --> returns the solar azimuth
--altitude() -->solar altitude
--zenith()--> solar zenith

-- After running the script (or when calling update()), the program will retrieve the desired
-- shading device position and luminance power. This will be the latest win1, win2 win3
-- or lum1 lum2 lum3 desired. If these variables do not exist, the last one will be
-- assumed.

-- Example 1: Leaving the window 1 on its 6th position of shading device, and turning off the
-- first luminaire set:

--			win1=6
--			lum1=0

-- Example 2: Turn the lights on between 7am and 8pm, and off any other time. Do not control the shades.

--			if hour() > 7 and hour()<20 then
--				lum1=1
--			else
--				lum1=0
--			end

-- Example 3: If irradiance in the south-facing facade is more than 300 W/m2, close the windows.
-- Then, update the illuminance levels, and see if we should turn the lights on.

--			if irradiance(0,-1,0) > 300 then
--				win1=1
--			else
--				win1=3
--			end
--			
--			update()
--			
--			if sensor(1) < 300 then
--				lum1=1
--			elseif sensor(1) > 1000 then
--				lum1=0
--			end
			

-- These strategies can be combined in order to create some complex control algorithms, and
-- even behavioral models.
