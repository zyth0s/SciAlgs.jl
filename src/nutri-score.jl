
# NUTRI-SCORE as developed by Sant'e France 2017
# Improved: 
# Development of a new front-of-pack nutrition label in France: the five-colour Nutri-Score
# Chantal Julia, Serge Hercberg

using TerminalMenus

# values /100 g ( or 100ml?)

food_types = ["water", "drink", "fats", "cheese", "other"]
food_type_menu = RadioMenu(food_types,pagesize=5)
type             = request("Choose type of food:",food_type_menu)
type = food_types[type]
print("Energy             [kJ/100 g]: "); energy           = parse(Float64,readline(stdin))
print("Fat                [ g/100 g]: "); fat              = parse(Float64,readline(stdin))
print("Saturated Fat      [ g/100 g]: "); sat_fat          = parse(Float64,readline(stdin))
print("Protein            [ g/100 g]: "); protein          = parse(Float64,readline(stdin))
print("Salt               [ g/100 g]: "); salt             = parse(Float64,readline(stdin))
print("Fruit & vegetables [ g/100 g]: "); fruit_vegetables = parse(Float64,readline(stdin))
print("Fibre              [ g/100 g]: "); fibre            = parse(Float64,readline(stdin))
print("Sugar              [ g/100 g]: "); sugar            = parse(Float64,readline(stdin))

# Example: Olive oil
#type             = "fats"
#energy           = 3381 # kJ
#fat              = 91.4    
#sat_fat          = 12.8
#protein          =  0.0
#salt             =  0.0
#fruit_vegetables =  0.0
#fibre            =  0.0
#sugar            =  0.0
# Result: 11 (D)

function energy_score(energy,type)
   if type == "drink"
      if     energy <=   0; pts=0
      elseif energy <=  30; pts=1
      elseif energy <=  60; pts=2
      elseif energy <=  90; pts=3
      elseif energy <= 120; pts=4
      elseif energy <= 150; pts=5
      elseif energy <= 180; pts=6
      elseif energy <= 210; pts=7
      elseif energy <= 240; pts=8
      elseif energy <= 270; pts=9
      elseif energy  > 270; pts=10
      end
   else
      if     energy <= 335;  pts=0
      elseif energy <= 670;  pts=1
      elseif energy <= 1005; pts=2
      elseif energy <= 1340; pts=3
      elseif energy <= 1675; pts=4
      elseif energy <= 2010; pts=5
      elseif energy <= 2345; pts=6
      elseif energy <= 2680; pts=7
      elseif energy <= 3015; pts=8
      elseif energy <= 3350; pts=9
      elseif energy >  3350; pts=10
      end
   end
   pts
end

function fat_score(fat,sat_fat)
   if fat > 100
      fat = 100
   end
   if sat_fat > 100
      sat_fat = 100
   end
   if type == "fats"
      ratio = 0
      if fat != 0
         ratio = round(sat_fat/fat*1e4)/1e2
      end
      if     ratio <  10; pts=0 
      elseif ratio <  16; pts=1 
      elseif ratio <  22; pts=2 
      elseif ratio <  28; pts=3 
      elseif ratio <  34; pts=4 
      elseif ratio <  40; pts=5 
      elseif ratio <  46; pts=6 
      elseif ratio <  52; pts=7 
      elseif ratio <  58; pts=8 
      elseif ratio <  64; pts=9 
      elseif ratio >= 64; pts=10 
      end
   else
      if     sat_fat <=  1; pts=0
      elseif sat_fat <=  2; pts=1
      elseif sat_fat <=  3; pts=2
      elseif sat_fat <=  4; pts=3
      elseif sat_fat <=  5; pts=4
      elseif sat_fat <=  6; pts=5
      elseif sat_fat <=  7; pts=6
      elseif sat_fat <=  8; pts=7
      elseif sat_fat <=  9; pts=8
      elseif sat_fat <= 10; pts=9
      elseif sat_fat >  10; pts=10
      end
   end
   pts
end

function sugar_score(sugar,type)
   if sugar > 100
      sugar = 100
   end
   if type == "drink"
      if     sugar <=  0  ; pts=0
      elseif sugar <=  1.5; pts=1
      elseif sugar <=  3  ; pts=2
      elseif sugar <=  4.5; pts=3
      elseif sugar <=  6  ; pts=4
      elseif sugar <=  7.5; pts=5
      elseif sugar <=  9  ; pts=6
      elseif sugar <= 10.5; pts=7
      elseif sugar <= 12  ; pts=8
      elseif sugar <= 13.5; pts=9
      elseif sugar >  13.5; pts=10
      end
   else
      if     sugar <=  4.5; pts=0
      elseif sugar <=  9  ; pts=1
      elseif sugar <= 13.5; pts=2
      elseif sugar <= 18  ; pts=3
      elseif sugar <= 22.5; pts=4
      elseif sugar <= 27  ; pts=5
      elseif sugar <= 31  ; pts=6
      elseif sugar <= 36  ; pts=7
      elseif sugar <= 40  ; pts=8
      elseif sugar <= 45  ; pts=9
      elseif sugar >  45  ; pts=10
      end
   end
   pts
end

function fibre_score(fibre)
   if fibre > 100
      fibre = 100
   end

   if     fibre <= 0.9; pts=0
   elseif fibre <= 1.9; pts=1
   elseif fibre <= 2.8; pts=2
   elseif fibre <= 3.7; pts=3
   elseif fibre <= 4.7; pts=4
   elseif fibre >  4.7; pts=5
   end
   pts
end

function protein_score(protein)
   if protein > 100
      protein = 100
   end
   if protein < 0
      protein = 0
   end
   if     protein <= 1.6; pts=0
   elseif protein <= 3.2; pts=1
   elseif protein <= 4.8; pts=2
   elseif protein <= 6.4; pts=3
   elseif protein <= 8  ; pts=4
   elseif protein >  8  ; pts=5
   end
   pts
end

function sodium_score(salt)
   if salt > 100
      salt = 100
   end
   sodium = salt*400
   if     sodium <=  90; pts=0
   elseif sodium <= 180; pts=1
   elseif sodium <= 270; pts=2
   elseif sodium <= 360; pts=3
   elseif sodium <= 450; pts=4
   elseif sodium <= 540; pts=5
   elseif sodium <= 630; pts=6
   elseif sodium <= 720; pts=7
   elseif sodium <= 810; pts=8
   elseif sodium <= 900; pts=9
   elseif sodium >  900; pts=10
   end
   pts
end

function fruit_score(fruit,type)
   if fruit > 100
      fruit = 100
   end
   if type == "drink"
      if     fruit <= 40; pts=0
      elseif fruit <= 60; pts=2
      elseif fruit <= 80; pts=4
      elseif fruit >  80; pts=10
      end
   else
      if     fruit <= 40; pts=0
      elseif fruit <= 60; pts=1
      elseif fruit <= 80; pts=2
      elseif fruit >  80; pts=5
      end
   end
   pts
end

function nutriscore_label(score,type)
   if type == "drink"
      if     score <= 1; score_label="B"
      elseif score <= 5; score_label="C"
      elseif score <= 9; score_label="D"
      else             ; score_label="E"
      end
   else
      if     score <= -1; score_label="A"
      elseif score <=  2; score_label="B"
      elseif score <= 10; score_label="C"
      elseif score <= 18; score_label="D"
      else              ; score_label="E"
      end
   end
   score_label
end

function nutri_score(type, energy, fat, sat_fat, protein,
                     salt, fruit_vegetables, fibre,
                     sugar) 
   score = 999
   nutri_score_label = "-"
   if type == "water"
      score_label = "A"
      return score, score_label
   end
   A = energy_score(energy,type)   +
       fat_score(fat,sat_fat)      + 
       sugar_score(sugar,type)     +
       sodium_score(salt)
   fruit_pts = fruit_score(fruit_vegetables,type)
   fibre_pts = fibre_score(fibre)
   C = fibre_pts     +
       protein_score(protein) +
       fruit_pts
    if (A >= 11) & (type != "cheese") & (fruit_pts < 5)
      score = A - fibre_pts - fruit_pts
   else
      score = A - C
   end
   score, nutriscore_label(score,type)
end

score, label = nutri_score(type, energy, fat, sat_fat, protein, salt, fruit_vegetables, fibre, sugar)
println("")
println("Result")
println("======")
println("The NUTRI-SCORE is $score with label $label")

