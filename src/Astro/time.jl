
# Box 2.2 in Fundamental Astronomy 6ed. by Karttunen et al.

function hms2h(hour,min,sec)
   hour + min/60 + sec/3600 # hours
end

function dms2d(deg,min,sec)
   deg + min/60 + sec/3600 # deg
end

function hour2deg(hour)
   hour*15.0
end

function deg2hour(deg)
   deg/15.0
end

function julian_date(y::Int,m::Int,d::Int)
   # returns the Julian date since noon (12:00 UT) w.r.t. 4713 B.C.
   m ∈ 1:12 || error("Wrong month m = $m.")
   d ∈ 1:31 || error("Wrong day d = $d.")
   367y  - ( 7(y + (m+9)÷12))÷4 -
   (3((y+(m-9)÷7)÷100+1)) ÷ 4 +
   275m÷9 + d + 1721029
end

function julian_date_Montenbruck1989(Y::Int,M::Int,D::Int,UT)
   # Montenbruck, O. (1989). Practical Ephemeris Calculations (Springer-Verlag, Heidelberg).
   # Described in Solar System Dynamics Appendix A.3
   # Y: year, M: month, D: day, UT: time of the day in hours wrt UT
   y = NaN; m = NaN; b = NaN
   if M < 3
      y = Y - 1
      m = M + 12
   else
      y = Y
      m = M
   end
   if Y < 1582
      B = -2
   end
   if Y == 1582 && 
      if M < 10
         B = -2
      end
      if M == 10 && D < 5
         B = -2
      end
      if M == 10 && D > 14
         B = floor(y/400) - floor(y/100)
      end
      if M > 10
         B = floor(y/400) - floor(y/100)
      end
   end
   if Y > 1582
      B = floor(y/400) - floor(y/100)
   end
   if y > 0
      return floor(365.25y) + floor(30.6001(m+1)) + B + 1720996.5 + D + UT/24
   else
      return trunc(365.25y-0.75) + floor(30.6001(m+1)) + B + 1720996.5 + D + UT/24
   end
end

function julian_date_inv_Montenbruck1989(j)
   # Montenbruck, O. (1989). Practical Ephemeris Calculations (Springer-Verlag, Heidelberg).
   frac(x) = modf(x)[1] # Fractional part
   a = floor(j+0.5)
   c = NaN
   if a < 2299161
      c = a + 1524
   else
      b = floor((a-1867216.25)/36524.25)
      c = a + b - floor(b/4) + 1525
   end
   d = floor((c-122.1)/365.25)
   e = floor(365.25d)
   f = floor((c-e)/30.6001)
   D = c - e - floor(30.6001f) + frac(j+0.5)
   M = f - 1 - 12floor(f/14)
   Y = d - 4715 - floor((7+M)/10)
   Int(Y),Int(M),D
end


function julian_date_J2000(y::Int,m::Int,d::Int)
   # returns the Julian date w.r.t. J2000.0 epoch
   julian_date(y,m,d) - 2451545.0
end

function julian_date_inv(J::Int)
   # J is the Julian date at noon (12:00 UT) w.r.t. 4713 B.C.
   a = J + 68569
   b = (4a)÷146097
   c = a - (146097b + 3)÷4
   d = (4000(c+1))÷1461001
   e = c - (1461d)÷4 + 31
   f = (80e)÷2447
   day = e - (2447f)÷80
   g = f÷11
   month = f + 2 - 12g
   year = 100(b-49) + d + g
   (year,month,day)
end

function weekday_from_julian(J::Int)
   if mod(J,7) == 0
      return "Monday"
   elseif mod(J,7) == 1
      return "Tuesday"
   elseif mod(J,7) == 2
      return "Wednesday"
   elseif mod(J,7) == 3
      return "Thursday"
   elseif mod(J,7) == 4
      return "Friday"
   elseif mod(J,7) == 5
      return "Saturday"
   elseif mod(J,7) == 6
      return "Sunday"
   end
end

function julian_century(J)
   # Time elapsed since January 1, 2000, in Julian centuries.
   # 2.48 in Fundamental Astronomy 6ed. by Karttunen
   (J - 2451545.0)/36525
end

function days_since_J2000(J)
   # J is the Julian date at noon (12:00 UT) w.r.t. 4713 B.C.
   J - 2451545.0
end

function GMST(T)
   # T is the Julian century
   24110.54841 + 8640184.812866T + 0.093104T^2 - 0.0000062T^3
end

function terrestial_time(TAI)
   # TAI is International Atomic Time [s]
   TAI + 32.184
end

testday = (1990,1,1)
julian_date_inv(julian_date(testday...)) == testday || error("Error")
weekday_from_julian(julian_date(testday...)) == "Monday" || error("Error")

#1946 Feb 4 10ʰ24ᵐ UT ≡ 1946 Feb 4 10.4ʰ UT
isapprox(julian_date_Montenbruck1989(1946,2,4,10.4), 2431855.933,atol=1e-3) || error("Error")
# 1954 June 10.25 UT ≡ 1954 June 10 6ʰ UT
julian_date_inv_Montenbruck1989(2434903.75) == (1954,6,10.25) || error("Error")



