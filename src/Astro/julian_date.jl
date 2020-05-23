
# Box 2.2 in Fundamental Astronomy 6ed. by Karttunen et al.

function julian_date(y::Int,m::Int,d::Int)
   367y  - ( 7(y + (m+9)÷12))÷4 -
   (3((y+(m-9)÷7)÷100+1)) ÷ 4 +
   275m÷9 + d + 1721029
end

function julian_date_inv(J::Int)
   # J is the Julian date at noon
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

testday = (1990,1,1)
julian_date_inv(julian_date(testday...)) == testday || error("Error")
weekday_from_julian(julian_date(testday...)) == "Monday" || error("Error")

