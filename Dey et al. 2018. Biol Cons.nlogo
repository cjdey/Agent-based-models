;This is a model of how polar bears and common eider ducks might respond to declines in Arctic sea ice.

;For model details see the full ODD model description associated with this model at github.com/cjdey/

;To run this model:
;-make sure you have the rngs extension downloaded (https://github.com/cstaelin/RNGS-Extension)
;-use Netlogo 5.3.1

;Then,
;Set desired model conditions -bears on? -bear decline? -clutch-size-climate-link? -breeding-propensity-climate-link?
;Set Up Landscape (may take a bit)
;Set Up Agents
;Go

;;The model is documented in

;Dey et al. 2017 Global Change Biology 23:1821:1831
;Dey et al. 2018 Biological Conservation 0:0


extensions [ rngs  profiler]

globals [

;;landscape related globals
island-seeds
number-islands
island-patches
day-of-year
hour-of-day
year
island-list

ice-breakup-date ;;; changes each year
historical-breakup-dates ;;; for 1990-2014
ice-freezeup-date ;;changes each year
historical-freezeup-dates ;;for 1989-2014
incidence-abundance
eb
first-mean-nesting-date
first-ice-breakup-date

;;bear related globals
bear-search-efficiency ;; a constant for Holling's disc equation.
bear-nest-handling-time ;; a constant for Holling's disc equation.
survival-rates  ;;a vector of annual bear survival rates
as  ;; the stable age-sex-distribution for bear agents
memN ;;memory coefficient for energy value of nests
memI ;; memory coefficient for islands
Qcrit-females
Qcrit-males

;;eider related globals
initial-number-nests
nest-clumpiness
beta ;;coefficient for variation in habitat quality
mean-nesting-date ;; changes each year
historical-mean-nesting-date  ;;for the start of the simulation
nesting-phenology-adjustment
annual-nest-failure-rate
initial-nest-distribution ;;; by island, to see changes from in future years
nest-distribution-by-cell
hatched-nests
active-nests
failed-nests
depredated-nests
annual-hatchling-survival-rate
annual-breeding-threshold
eider-senescence-rate
initial-eider-age-dist
eider-senescence-a
eider-senescence-b
eider-adult-survival-rate

;;; model outputs
depredated-nests-by-year
adults-by-year
juveniles-by-year
fledglings-by-year
ice-breakup-by-year
ice-freezeup-by-year
year-list
mean-nesting-date-by-year
clutch-size-by-year
var-nests-per-cell
mean-nest-ycor
var-nest-ycor
nests-in-large-colonies
nests-per-island
nests-per-island-by-year
breeding-propensity-by-year
hatchling-survival-by-year
mean-age-adults-by-year
mean-age-adults-and-juveniles-by-year


;; globals for observations
island-size
island-visits
island-distance
island-nests-here

bear-days-by-year
bear-days
bear-days-tick
]

breed [nests nest]
breed [bears bear]

bears-own [
  sex
  age-class
  age
  daily-E
  body-mass
  spring-body-mass
  predicted-body-mass ;;at freezeup date
  body-length
  onshore-date  ;;;date they appear on landscape each year
  energetic-stress
  state                ;;;what behavior they are going to do
  expected-E-per-nest
  potential-nest-intake
  I-should-local-forage?
  attendant-cub?
  attendant-yearling?
  where-I-sense-nests
  memory-islands
  memory-nests
  memory-islands-this-year
  memory-nests-this-year
  potential-targets
  target ;;;where I am heading, if swimming long distance
  I-have-visited-best-island?
  I-should-rest?
  hours-awake
  hours-resting
  dmax ;;max distance I can move in a time-step, determined by sex and parental status
  mass-loss-final
  mass-loss-pred
       ]

nests-own [
  incubate-time-remaining
  last-year-fate
  last-year-location
  nest-energy
  clutch-size
  nest-initiation-date
  nest-failure-rate
  age
  age-class
  breeding-this-year?
  qual
  ]

patches-own [
  island-ID
  habitat-quality
  bear-visits
  incidence-abundance-grid
    ]



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;SETUP LANDSCAPE;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-landscape
clear-all
ask patches [ set pcolor blue ]
rngs:init rngs:set-seed 1 232  ;; set up a random-seed

;;;create some islands
set island-seeds 500
let mylist n-values island-seeds [ ? + 1 ]
let yrange max-pycor - min-pycor - 61   ;;; the possible range for island centers
  foreach mylist [
        ask patch random-xcor ( (max-pycor - 61) - (yrange * ( rngs:rnd-beta 1 0.63 1.927) ) )
         [set pcolor green set island-ID ?
          let island-radius round sqrt ( (random-gamma 0.22243 0.003085) / 3.14 )
          let nearby-island-patches [island-ID] of other patches in-radius (island-radius + 1)  with [pcolor = green]
          ifelse (length nearby-island-patches > 0)
             [let new-ID one-of nearby-island-patches
             ask patches with [member? island-ID nearby-island-patches] [set island-ID new-ID]
             ask patches in-radius (island-radius) [set pcolor green set island-ID new-ID ] ]
             [ask patches in-radius (island-radius) [set pcolor green set island-ID [island-ID] of myself]]
  ]]

ask (patches with [pycor >= max-pycor - 45] )  [set pcolor grey set island-ID 0]   ;; this area will be the mainland, so remove any island areas that extend into mainland
ask (patches with [pycor = max-pycor - 46] )  [set pcolor blue set island-ID 0]   ;; make sure there is at least a one cell gap between any islands and the mainland ;;

set island-patches patches with [pcolor = green]

;renumber the islands from 1 to number-islands
let a sort remove-duplicates [island-ID] of island-patches
set number-islands length a
let b n-values number-islands [ ? + 1 ]
  foreach a [ask island-patches with [ island-ID = ? ] [ set island-ID item (position ? a) b] ]

set eb island-patches with [island-ID = 1]

;;;;;;set permanent climatic conditions
set historical-breakup-dates [190 182 194 180 177 169 182 173 166 169 169 168 180 171 178 170 158 175 171 176 160 172 164 173 173]
set historical-freezeup-dates [295 288 296 302 304 310 311 312 298 314 303 308 314 319 320 308 313 327 311 307 304 335 315 314 309 308]  ;;; one additional freezeup date compared to breakup dates (for year before first year)
set historical-mean-nesting-date  165 ;;; This is the date of first nesting (not mean nesting initiation date)
set nesting-phenology-adjustment 0.2

;;;setup model output recorders
set depredated-nests-by-year [] set year-list [] set ice-breakup-by-year [] set ice-freezeup-by-year [] set mean-nesting-date-by-year []
set mean-nest-ycor [] set var-nest-ycor [] set nests-in-large-colonies [] set var-nests-per-cell [] set adults-by-year [] set juveniles-by-year [] set fledglings-by-year [] set nests-per-island [] set nests-per-island-by-year []
set breeding-propensity-by-year [] set hatchling-survival-by-year [] set incidence-abundance [] set mean-age-adults-by-year [] set mean-age-adults-and-juveniles-by-year [] set clutch-size-by-year []

set island-list n-values number-islands [? + 1]
reset-ticks
end


 ;;;;;;;;;;;;;;;;;;;;;;
 ;;;;;;; SETUP AGENTS ;;;;;;;;
 ;;;;;;;;;;;;;;;;;;;;;

to setup-agents
  clear-turtles
   setup-nests
   if bears-on? = TRUE [  setup-bears ]
   set hour-of-day 0   set year 1   set day-of-year 121
end

to setup-nests
  set-default-shape nests "bird"
  set initial-number-nests 5000 set nest-clumpiness 0.40 set beta 0.2 set memI 0.5 set memN 0.5 set annual-nest-failure-rate 0.15 set hatched-nests 0
  set eider-adult-survival-rate 0.9225
  set eider-senescence-a 4.74E-11 set eider-senescence-b 7.431
calculate-eider-age-distribution

  rngs:init rngs:set-seed 1 77

let p ( (100 * nest-clumpiness) / number-islands)
        create-nests  initial-number-nests [
           set color brown
           set size 2  ;; easier to see
           set age one-of initial-eider-age-dist
           set qual random-normal 0 1
           ifelse age >= 3 [set age-class "adult"] [set age-class "juvenile"]
           ifelse random-float 1 > annual-nest-failure-rate [set last-year-fate 1] [set last-year-fate 0]
          hide-turtle
          set incubate-time-remaining 0
          set mean-nesting-date 0
          set last-year-location ( (rngs:rnd-negbinomial 1  1  (1 - p)) + 1 )      ;;; need to use prob of failure (not success)
  ]

 ask nests [
  let a last-year-location
   let where-to-nest island-patches with [ island-ID = a  and count turtles-here < 250]  ;;; 250 is a max density per cell (equivalent to 1000 nests / ha)
  ifelse any? where-to-nest
  [move-to one-of where-to-nest ]
  [ move-to one-of island-patches with [count turtles-here < 250 ] set last-year-location [island-ID] of patch-here]
         ]

set initial-nest-distribution [island-ID] of nests


;;;setup habitat quality ie nest success rate of cells
let mean-island mean [island-ID] of nests
 ask island-patches [
     let x island-ID
     let y length filter [ ? = x ] initial-nest-distribution
     set habitat-quality (annual-nest-failure-rate -   (beta / 100) * (mean-island - x) ) ;;beta is decrease in nest failure rate per island rank (in percentage points)
     if habitat-quality > 1 [set habitat-quality 0.99]
   ]

end

 to setup-bears
  set bear-search-efficiency 1  ;; equivalent to the number of patches per tick a bear can fully search. needed for the Holling disc eqn
  set bear-nest-handling-time 0.01666667 ;; equal to 1/60th of a tick (30s per nest)
  set Qcrit-males 51.1
  set Qcrit-females 38.6
  ;;;;calculate the stable bear age and sex distribution ;;;need to do this to define initial sex and ages of bear agents
   set as[]
   set survival-rates [0.81 0.90 0.75 0.86 0.93 0.81 0.595]   ;; male subadults, adult, old adults then female subadult, adult, old adults, then cubs of both sexes
   calculate-stable-age-sex-distribution

  set-default-shape bears "bear"
   create-bears 40
  [
     let r random length as
     set age item 0 item r as ;;; sample from the stable age/sex distribution
     set sex item 1 item r as
 ifelse age < 5 [set age-class "subadults"] [set age-class "adult"]
    set color white set size 5
 let i sort remove-duplicates ([island-ID] of island-patches)
    set memory-islands i ;;; a list of islands
    set memory-nests n-values (length memory-islands) [0] ;;;expected nests for each island
    set memory-islands-this-year [ ]
    set memory-nests-this-year [ ]
    set target nobody
    set I-have-visited-best-island? FALSE
    set expected-E-per-nest ((random-float 1.61) + 2.32) ;; random value between the minimum E content per nest (old nest) and the max (fresh nests)

    set attendant-cub? FALSE set attendant-yearling? FALSE
    if age-class ="adult" and sex ="female" [ifelse random-float 1 < 0.33 [set attendant-cub? TRUE] [if random-float 1 > 0.5 [set attendant-yearling? TRUE]]]
        hide-turtle
  ]

 update-plots
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;; GO  ;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to go

 if day-of-year = 121 and hour-of-day = 0 [new-year-procedures]     ;;;;do some stuff at the start of each year
 if count nests with [breeding-this-year? = TRUE] = 0 [ stop]
if bears-on? = TRUE [
                  ask bears [
            if hour-of-day = 0 [ bears-appear ] ;;occurs once per year (at the bear onshore date)
            if hidden? = FALSE [                ;; if they are active
              check-energetic-stress
              update-state
              behave]
            if hidden? = FALSE and hour-of-day = 0 [ convert-daily-E-to-mass]
   ]
        ;record-bear-days
        ;if year = 21 [  record-bear-visits ]
            ]

if hour-of-day = 0 [ ask nests [
        if hidden? = FALSE [            ;;if they are active
              random-fail
              nest-success
              incubate]
        if day-of-year = nest-initiation-date and breeding-this-year? = TRUE  [     initiate-nest    ]
  ]]

count-time
  if year > 50 [ record-incidence-abundance stop]
  tick
end

;;;;;;;;;;;;;;;;;;;
;;;; PATCH/LANDSCAPE PROCEDURES;;;
;;;;;;;;;;;;;;;;;;;
to new-year-procedures

    setup-climate
    clear-drawing
    set hatched-nests 0 set failed-nests 0 set depredated-nests 0 set active-nests 0
    if year = 1 [set first-mean-nesting-date mean-nesting-date set first-ice-breakup-date ice-breakup-date]
    calculate-hatchling-survival-rate
    calculate-annual-breeding-threshold
        if bears-on? = TRUE [
    ask bears
     [ hide-turtle
        age-and-death
        age-my-memory
                ]

if bear-decline? = TRUE [check-bear-death]

     ask bears [
       move-to-ice-edge
       calculate-arrival-date
       calculate-spring-condition
       determine-attendant-offspring-and-dmax
       set daily-E 0
       set hours-awake random 48 / 2]
    ]

     ask nests [
          if year > 1 [get-older check-death]
        check-breeding-this-year
            if breeding-this-year? = TRUE
         [calculate-initiation-date
           calculate-clutch-size
           if year > 1 [relocate-nest]]
                    ]

    set nest-distribution-by-cell [count nests-here] of island-patches
    ask nests [set-nest-failure-rate ]

   ifelse count nests with [breeding-this-year? = TRUE] = 0 [set day-of-year 100] [ set day-of-year min [nest-initiation-date] of nests with [breeding-this-year? = TRUE] ]


end

to setup-climate
  ifelse year < 25 [set ice-breakup-date item (year - 1) historical-breakup-dates]
           [set historical-breakup-dates lput ice-breakup-date historical-breakup-dates set ice-breakup-date round (181.86 - (0.635 * year) + random-normal 0 6.994) ]
  ifelse year < 25 [set ice-freezeup-date item (year) historical-freezeup-dates]
           [set historical-freezeup-dates lput ice-freezeup-date historical-freezeup-dates set ice-freezeup-date round (300.33 + (0.727 * year) + random-normal 0 8.31) ]
    let spring-advance (mean historical-breakup-dates) - (ice-breakup-date) ;; how much is the spring advanced from a typical year
  set mean-nesting-date round (historical-mean-nesting-date - (spring-advance * nesting-phenology-adjustment) )
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;            NEST PROCEDURES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calculate-initiation-date
  set nest-initiation-date round (random-normal (mean-nesting-date) 7.14 )
end

to calculate-clutch-size
  let intercept (4.12 - (-0.027009 * first-mean-nesting-date))
  let x 0
  ifelse clutch-size-climate-link? = TRUE  [ set x (intercept + (-0.027009 * mean-nesting-date)) ] [ set x 4.12] ;;the expected value given the nesting date
   set clutch-size (rngs:rnd-poisson 1 (x - 1 + random-normal 0 0.3039346)  ) + 1
end

to calculate-hatchling-survival-rate
 set annual-hatchling-survival-rate (rngs:rnd-beta 1 1.272 5.423 )
 end


to calculate-annual-breeding-threshold
ifelse breeding-propensity-climate-link? = TRUE
[ set annual-breeding-threshold (inverse-logit (1.65391 + (-0.01111 * ice-breakup-date) + random-normal 0 0.4478635))]
[ set annual-breeding-threshold (inverse-logit (1.65391 + (-0.01111 * first-ice-breakup-date) + random-normal 0 0.4478635))]   ;;;conditions during year 1

end


 to check-breeding-this-year
  ifelse age-class = "adult"
  [ifelse random-float 1 < annual-breeding-threshold [set breeding-this-year? TRUE][set breeding-this-year? FALSE set clutch-size 0]]
   [set breeding-this-year? FALSE set clutch-size 0]
end

to check-death
   ifelse age-class = "hatchling"
   [ if random-float 1 > eider-adult-survival-rate ^ (11 / 12) [die] ]  ;;juvenile survival rate is same as adults but they have already survived 1 month (hatch -> 30 days)
   [
   ifelse age <= 2 [if random-float 1 > eider-adult-survival-rate [ die]]
                    [if random-float 1 < ( (1 - eider-adult-survival-rate) + (eider-senescence-a * ((age - 2) ^ eider-senescence-b))) [die]]
     ]
 end

to get-older
  set age age + 1
  if age < 3 [set age-class "juvenile"]
  if age >= 3 [set age-class "adult"]
end

to relocate-nest
  let dispersal-distance 0
       ifelse last-year-fate = 0
       [set dispersal-distance round ((rngs:rnd-lognormN 1 (3.72 ) 1.21) / 50 )   ]
       [set dispersal-distance round ((rngs:rnd-lognormN 1 (2.80 ) 1.40) / 50 )  ]
 if dispersal-distance > 1 [
           let y patches in-radius dispersal-distance
           let new-nest-spot (y with [pcolor = green and count nests-here < 250 ])
           if any? new-nest-spot  [move-to max-one-of new-nest-spot [distance myself]
       ]]
 end

to set-nest-failure-rate
     set nest-failure-rate habitat-quality
end

to initiate-nest
    show-turtle
    set active-nests active-nests + 1
    set incubate-time-remaining 26       ;;;; 26 days of incubation
    set nest-energy clutch-size * 0.984 ;;;;energy in a fresh egg is 1.087 MJ * digestibility of .905
end

to random-fail
 if random-float 100 > (100 * ( (1 - nest-failure-rate) ^ ( 1 / 26))) [
    hide-turtle
    set incubate-time-remaining 0
    set last-year-fate 0
    set nest-energy 0
    set failed-nests failed-nests + 1
    set active-nests active-nests - 1]
end

to nest-success
  if incubate-time-remaining = 1 [
    hide-turtle
    set incubate-time-remaining 0
    set hatched-nests hatched-nests + 1
    set active-nests active-nests - 1
    set nest-energy 0
    set last-year-fate 1
    create-ducklings
    ]
end

to create-ducklings
  let number-female-hatchlings random-binomial clutch-size 0.5 ;; assumes all eggs hatch i.e. no partial nest predation
  let number-female-hatchlings-survived (random-binomial number-female-hatchlings annual-hatchling-survival-rate)
  hatch-nests number-female-hatchlings-survived ;;;create new eiders
  [set age 0 set age-class "hatchling"
  set last-year-fate 0
  set hidden? TRUE
  set qual random-normal 0 1
  set clutch-size 0
  set breeding-this-year? FALSE
  ]
end

to incubate   if hidden? = FALSE [
  set incubate-time-remaining incubate-time-remaining - 1
  ifelse incubate-time-remaining > (26 / 3)
            [set nest-energy 0.984 * clutch-size]   ;;nest energy is the same for the first 2/3 of incubation
            [set nest-energy ((incubate-time-remaining * 0.050416419) + 0.5804) * clutch-size]  ;;then decreases linearly over the last 1/3
   ]
end

to calculate-eider-age-distribution
set initial-eider-age-dist n-values 10000 [(random 20) + 1]
 repeat 1000 [  set initial-eider-age-dist map
      [ifelse-value (? <= 2) [ifelse-value (random-float 1 > eider-adult-survival-rate) [ 1 ] [(? + 1)]]
                      [ifelse-value (random-float 1 < ( (1 - eider-adult-survival-rate) + (eider-senescence-a * ((? - 2) ^ eider-senescence-b))))     [ 1 ] [(? + 1)]]

]                      initial-eider-age-dist   ]
  end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;BEAR PROCEDURES;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to age-and-death
ifelse (sex = "male")
;;if they are males
 [ifelse  (age < 5)   [ifelse (random-float 1 > item 0 survival-rates) [reset-bear-agent] [set age age + 1 ] ]
      [ ifelse (age <= 19) [ifelse (random-float 1 > item 1 survival-rates) [reset-bear-agent] [set age age + 1 ] ]
        [ifelse (random-float 1 > item 2 survival-rates ) [reset-bear-agent] [set age age + 1 ] ] ]
  ]
 ;;if they are females
  [ifelse (age < 5)   [ifelse (random-float 1 > item 3 survival-rates) [reset-bear-agent] [set age age + 1 ] ]
      [ ifelse (age <= 19) [ifelse (random-float 1 > item 4 survival-rates) [reset-bear-agent] [set age age + 1 ] ]
        [ifelse (random-float 1 > item 5 survival-rates ) [reset-bear-agent] [set age age + 1 ] ] ]
    ]

    ifelse age >= 5 [set age-class "adult"] [set age-class "subadult"]
end

to check-bear-death
if (member? year [26 28 30 32 34 36 38 40 42 44 46 48]) [ask one-of bears [die]]
end

to reset-bear-agent
  set age 3 set age-class "subadult"
  set sex random-sex
  set memory-nests n-values (length memory-islands) [0]
  set memory-islands-this-year [ ]       set memory-nests-this-year [ ]
  set I-have-visited-best-island? FALSE
  set attendant-cub? FALSE set attendant-yearling? FALSE
  set expected-E-per-nest ((random-float 1.61) + 2.32)
end

to move-to-ice-edge
  setxy random-xcor (min-pycor + random 50)
  set heading 0
end

to calculate-arrival-date
set onshore-date round (ice-breakup-date + random-normal 24.6 5.5)
end

to calculate-spring-condition
  ifelse sex = "male"   ;;;;set a new body length
  [set body-length (237 * (1 - exp(1)^(-0.381 * (age + 1.219)))) / 100 ] ;;length from VonB curve
  [set body-length (198 * (1 - exp(1)^(-0.530 * (age + 1.166)))) / 100 ] ;;length from VonB curve

  ifelse sex = "male"
  [  set body-mass (579 * (1 - exp(1)^(- 0.206 * (age - (-3.383)))) ^ 3) ;;vonB curve for males (for mass)
    set body-mass random-normal body-mass (67.6) ]

    [ set body-mass (255 * (1 - exp(1)^(- 0.605 * (age - (-0.944)))) ^ 3) ;; vonB curve for females (for mass)
    set body-mass random-normal body-mass (55.3)]
  if body-mass < 100 [calculate-spring-condition]

  set hours-awake random 16
  set hours-resting 0
  set spring-body-mass body-mass
end

to determine-attendant-offspring-and-dmax
  if sex = "male" [set dmax 12]
  if sex = "female" and age-class = "subadult" [ set dmax 20 ]
  if sex = "female" and age-class = "adult" [
    if attendant-cub? = TRUE [set attendant-cub? FALSE ifelse random-float 1 < item 2 survival-rates [set attendant-yearling? TRUE set dmax 18] [set attendant-yearling? FALSE set dmax 20]]
    if attendant-yearling? = TRUE [ set attendant-cub? FALSE set attendant-yearling? FALSE set dmax 20]
    if attendant-cub? = FALSE and attendant-yearling? = FALSE [set attendant-cub? TRUE set dmax 17]
  ]
end

to bears-appear
  if day-of-year = onshore-date [ show-turtle set daily-E 0 ]
end

to check-energetic-stress
    set predicted-body-mass (body-mass - ( (1 / 39.3) * energy-cost-resting * 48 * (item (year - 1) historical-freezeup-dates - day-of-year) ) )
    ifelse sex = "male"
         [ ifelse ( predicted-body-mass / body-length ^ 2) <= Qcrit-males ;
                                             [set energetic-stress "high"]  [set energetic-stress "low" ] ]
         [ ifelse ( predicted-body-mass / body-length ^ 2) <= Qcrit-females
                                             [set energetic-stress "high"]  [set energetic-stress "low"]]
end

to update-state
  check-rest ifelse I-should-rest? = TRUE [set state "rest" ]
      [check-local-forage ifelse I-should-local-forage? = TRUE [set state "forage"]
        [ifelse any? nests-I-can-smell [set state "move-to-nearby-nests" ]
          [ ifelse energetic-stress = "low" [set state "migrate-towards-mainland" ]
             [ifelse max memory-nests > 1 [check-working-memory-best-island ifelse I-have-visited-best-island? = TRUE [set state "search-for-food"] [set state "memory-informed-search"]]
              [set state "search-for-food"]]
      ]]]
end

to check-rest
  ifelse pcolor = blue
     [set I-should-rest? FALSE]
     [ifelse energetic-stress = "high"
                  [ ifelse  hours-awake > (24 * 0.266)  [set I-should-rest? TRUE ] [set I-should-rest? FALSE ] ]
                  [ ifelse  hours-awake > (24 * 0.145)  [set I-should-rest? TRUE ] [set I-should-rest? FALSE ] ]
     ]
end

to check-local-forage ;;check if it is worth it to forage here
  ifelse pcolor = green
 [ let local-active-nests count nests-here with [hidden? = FALSE]
  set potential-nest-intake round (  (bear-search-efficiency * local-active-nests) / (1 + (bear-search-efficiency * bear-nest-handling-time * local-active-nests))) ;number of nests I will get according to Hollings Disc
  ifelse (potential-nest-intake * (expected-E-per-nest) ) > energy-cost-active [set I-should-local-forage? TRUE] [set I-should-local-forage? FALSE];; if local foraging would provide a net E gain
  ]
 [set I-should-local-forage? FALSE]
end

to check-working-memory-best-island
   let x position (max memory-nests) memory-nests
    let y item x memory-islands
    ifelse member? y memory-islands-this-year [set I-have-visited-best-island? TRUE] [set I-have-visited-best-island? FALSE]
end

to behave
    if state = "rest" [rest set label "resting" ]
    if state = "forage" [forage set label "foraging" ]
    if state = "move-to-nearby-nests" [move-to-nearby-nests set label ""]
    if state = "migrate-towards-mainland" [migrate-towards-mainland set label ""]
    if state = "memory-informed-search" [memory-informed-search set label "memory-informed search"]
    if state = "search-for-food" [search-for-food set label "random-search"]
end

to rest
   set hours-resting hours-resting + 0.5
   ifelse energetic-stress = "high"
           [if hours-resting > (24 - (24 * 0.266)) [set hours-awake 0 set hours-resting 0]]
           [if hours-resting > (24 - (24 * 0.145)) [set hours-awake 0 set hours-resting 0]]
   set daily-E (daily-E - energy-cost-resting )   ;;;lose a little bit of energy while resting
end

to forage
     let prey n-of potential-nest-intake nests-here with [hidden? = FALSE]
     set daily-E (daily-E +  (sum [nest-energy] of prey)  - energy-cost-active )   ;;; bear gains the energy from nests but loses the energy cost of foraging
     set expected-E-per-nest (expected-E-per-nest * (1 - memN) + (memN * ((sum [nest-energy] of prey) / potential-nest-intake)))  ;; update expected-E-per-nest via linear operator model
      ask prey [
       hide-turtle
       set incubate-time-remaining 0
       set nest-energy 0
       set last-year-fate 0  ]
     set depredated-nests depredated-nests + potential-nest-intake
     set active-nests active-nests - potential-nest-intake
update-memory ;;remember there are nests here
set hours-awake hours-awake + 0.5
end

to move-to-nearby-nests
  update-memory  ;;remember there are no nests here
  let x min-one-of nests-I-can-smell [distance myself]
  face x
  move-to x
  set daily-E (daily-E - energy-cost-active ) ;; E cost of movement
  set hours-awake hours-awake + 0.5
end

to migrate-towards-mainland
   update-memory
    if pcolor = green [  set heading 0 right random-normal 0 45 walk set daily-E (daily-E - ( energy-cost-active )) set hours-awake hours-awake + 0.5 ]  ;;cost of movement
    if pcolor = blue  [ set heading 0 right random-normal 0 45 swim set daily-E (daily-E - ( energy-cost-active )) set hours-awake hours-awake + 0.5  ] ;;cost of movement
    if pcolor = grey [ set state "rest" rest ]  ;;; if I'm on the mainland rest
 end

to search-for-food
  update-memory  ;; remember there are no nests here
  ifelse pcolor = green
  [right random-normal 0 45 walk]
  [right random-normal 0 45 swim]
  set daily-E (daily-E - ( energy-cost-active ))
  set hours-awake hours-awake + 0.5
end

to memory-informed-search
  update-memory  ;;remember there are no nests here
  let x position (max memory-nests) memory-nests
  let y item x memory-islands  ;; y is the island-ID of the best island from memory

  ifelse island-ID = y [search-for-food] ;;;; If I am on the target island, search for food
                       [ifelse pcolor = green
                             [face one-of island-patches with [island-ID = y] right random-normal 0 45 walk]
                             [face one-of island-patches with [island-ID = y] right random-normal 0 45 swim]
                        ]

 set daily-E (daily-E - (energy-cost-active) )  ;;E cost of movement
 set hours-awake hours-awake + 0.5
end

to swim
       let x dmax
       ifelse target != nobody
          [ifelse distance target > x [face target move] [move-to target ] ]  ;;if I already have a target go towards it
          [ check-for-visible-islands ifelse any? potential-targets
             [ifelse any? potential-targets with [distance myself < x]  [move-to min-one-of potential-targets [distance myself] ] [set target min-one-of potential-targets [distance myself] face target move]] ;;if I can see an island that I haven't visited head towards it
             [ ifelse can-move? x and ( (heading < 90) or (heading > 270) or (ycor > 0) ) [forward x] [right 180 forward x] ] ;; if I can't see any islands I haven't visited, correlated random swim. Prevented from swimming out to sea
  ]
end

to check-for-visible-islands
  let y patches in-cone 60 180
  let x memory-islands-this-year
  set potential-targets y with [pcolor = green and not member? island-ID x]
end

to walk
  set potential-targets nobody   set target nobody  ;;;reset these since I am on land now
  ifelse can-move? dmax [] [right 180]  ;; if I can't move it (i.e. side of map, then turn around)
 move
 end

to move
   let move-step 0
 while [not any? nests-I-can-smell and move-step < dmax ] ;;;smell while moving, if I smell anything, then go there
 [ fd 1 set move-step (move-step + 1)]
  if any? nests-I-can-smell [
    let y min-one-of nests-I-can-smell [distance myself]
    face y
    move-to y]
end

to convert-daily-E-to-mass
      set body-mass (body-mass + (daily-E / 39.3) )
      set daily-E 0
end

to update-memory
 if pcolor = green [  ;if I'm on an island
   set target nobody set potential-targets nobody ;;erase my target
  ifelse state = "forage"[
   ifelse member? island-ID memory-islands-this-year   ;;;if this island is already in my memory
  [ set memory-nests-this-year
     replace-item (position island-ID memory-islands-this-year) memory-nests-this-year   ;;;replace the corresponding value in memory-nests-this-year
     (item (position island-ID memory-islands-this-year) memory-nests-this-year + potential-nest-intake) ]   ;; by adding the number of new nests I consumed on this island
  [ set memory-nests-this-year lput potential-nest-intake memory-nests-this-year ;; Otherwise make a new entry in each list ;;with the number of nests I ate
    set memory-islands-this-year lput island-ID memory-islands-this-year]
  ]
 [ifelse member? island-ID memory-islands-this-year   ;;;if this island is already in my memory
  [  ]   ;; do nothing (don't update the number of nests because I'm not eating any)
  [ set memory-nests-this-year lput 0 memory-nests-this-year ;; Otherwise make a new entry in each list ;; with 0 nests
    set memory-islands-this-year lput island-ID memory-islands-this-year]
 ]
 ]
end

to calculate-mass-loss
  ifelse onshore-date < day-of-year [set mass-loss-final ( (body-mass - spring-body-mass) / (day-of-year - onshore-date) )
   set mass-loss-pred ( (predicted-body-mass - spring-body-mass) / (ice-freezeup-date - onshore-date) ) ]
  [set mass-loss-final 9999 set mass-loss-pred 9999]
end


to age-my-memory
  set target nobody
  set potential-targets nobody
  set memory-nests map [? * (1 - memI) ]  memory-nests ;;;first degrade memory (age)
  (foreach memory-islands-this-year memory-nests-this-year
    [let x position ?1 memory-islands
       set memory-nests (replace-item x memory-nests (item x memory-nests + (?2 * memI) ))   ;then add any new information gained that year
      ] )
  set memory-islands-this-year [ ]  set memory-nests-this-year [ ]
  set I-have-visited-best-island? FALSE
end

to calculate-stable-age-sex-distribution ;;for bears
repeat 10000 [  let x list (random 25) random-sex set as lput x as ]
repeat 1000 [  set as map
  [   ifelse-value (item 1 ? = "male") [
   ifelse-value (item 0 ? < 5)   [ifelse-value (random-float 1 > item 0 survival-rates) [list 3 random-sex] [list (item 0 ? + 1) (item 1 ?)] ]         ;;;subadults
      [ ifelse-value (item 0 ? <= 19) [ifelse-value (random-float 1 > item 1 survival-rates) [list 3 random-sex] [list (item 0 ? + 1) (item 1 ?)] ]     ;;prime age adults
             [ifelse-value (random-float 1 > item 2 survival-rates ) [list 3 random-sex] [list (item 0 ? + 1) (item 1 ?)] ] ]  ;;;older adults
  ]
    [ifelse-value (item 0 ? < 5)   [ifelse-value (random-float 1 > item 3 survival-rates) [list 3 random-sex] [list (item 0 ? + 1) (item 1 ?)] ]
      [ ifelse-value (item 0 ? <= 19) [ifelse-value (random-float 1 > item 4 survival-rates) [list 3 random-sex] [list (item 0 ? + 1) (item 1 ?)] ]
             [ifelse-value (random-float 1 > item 5 survival-rates ) [list 3 random-sex] [list (item 0 ? + 1) (item 1 ?)] ] ]
    ]
  ] as  ]
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;           Time                       ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to count-time
   ifelse hour-of-day < 23.5
    [set hour-of-day hour-of-day + 0.5] [set hour-of-day 0]

  if hour-of-day = 0 [set day-of-year day-of-year + 1 ]

  if day-of-year > 170 and count nests with [hidden? = FALSE] < 1 [
    record-results ;;record the results from last season
    set year year + 1
    set hour-of-day 0
    set day-of-year 121
        ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;           Model Outputs                          ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to record-results  ;;end of year results
     set mean-age-adults-by-year lput (mean [age] of nests with [age-class = "adult"]) mean-age-adults-by-year
     set mean-age-adults-and-juveniles-by-year lput (mean [age] of nests with [age-class = "adult" or age-class = "juvenile"]) mean-age-adults-and-juveniles-by-year
     set fledglings-by-year lput (count nests with [age-class = "hatchling"]) fledglings-by-year
     set juveniles-by-year lput (count nests with [age-class = "juvenile"]) juveniles-by-year
     set adults-by-year lput (count nests with [age-class = "adult"]) adults-by-year
     set year-list lput year year-list


     set clutch-size-by-year lput ( mean [clutch-size] of nests with [clutch-size > 0]) clutch-size-by-year

     set depredated-nests-by-year lput depredated-nests depredated-nests-by-year
     set ice-breakup-by-year lput ice-breakup-date ice-breakup-by-year
     set ice-freezeup-by-year lput ice-freezeup-date ice-freezeup-by-year
     ifelse count nests = 0 [set mean-nesting-date-by-year lput "NA" mean-nesting-date-by-year] [set mean-nesting-date-by-year lput mean ([nest-initiation-date] of nests) mean-nesting-date-by-year]
     ifelse count nests = 0 [set mean-nest-ycor lput "NA" mean-nest-ycor] [set mean-nest-ycor lput mean ([ycor] of nests) mean-nest-ycor]


    ;set var-nest-ycor lput variance ([ycor] of nests) var-nest-ycor
    set nests-in-large-colonies lput (count nests with [island-ID < 11]) nests-in-large-colonies
    ifelse count nests = 0 [set var-nests-per-cell lput "NA" var-nests-per-cell] [set var-nests-per-cell lput (variance nest-distribution-by-cell) var-nests-per-cell]
    set breeding-propensity-by-year lput (count nests with [age-class = "adult" and breeding-this-year? = TRUE] / count nests with [age-class = "adult"]) breeding-propensity-by-year
    set hatchling-survival-by-year lput annual-hatchling-survival-rate hatchling-survival-by-year

    let x [ ]
    foreach island-list [set x lput count nests with [island-ID = ? and breeding-this-year? = TRUE] x]
    set nests-per-island-by-year lput x nests-per-island-by-year
end

to record-incidence-abundance
let cell-list [1 2 3 4 5 6 7 8 9 10]
foreach cell-list [
  let y ?
let x remove-duplicates [island-ID] of (island-patches with [pxcor >= (-1000 + ((y - 1) * 200)) and pxcor < (-800 + ((y - 1) * 200)) and pycor > -125 and pycor <= 75])
foreach x [  ask island-patches with [island-ID = ?] [set incidence-abundance-grid y]]
 foreach x [if any? (island-patches with [island-ID = ? and not (pxcor >= (-1000 + ((y - 1) * 200)) and pxcor < (-800 + ((y - 1) * 200)) and pycor > -125 and pycor <= 75)])
                                  [ask island-patches with [island-ID = ?] [set incidence-abundance-grid 0] ]]

let z remove-duplicates [island-ID] of island-patches with [incidence-abundance-grid = y] ;;;list of islands within that cell
let i [] foreach z [ifelse any? nests with [island-ID = ? and breeding-this-year? = TRUE] [set i lput 1 i][set i lput 0 i] ]
let a [] foreach z [set a lput (count nests with [island-ID = ? and breeding-this-year? = TRUE]) a ]
let l list ((sum i) / (length z))  ;;incidence
           ((sum a) / (length z))  ;;abundance

set incidence-abundance lput l incidence-abundance

]


end


;;;;;;;;;;;reporters

to-report energy-cost-resting
  let baseline-cost (body-mass * (0.002438 / 2) )
ifelse attendant-cub? = FALSE and attendant-yearling? = FALSE
[report baseline-cost]
[ifelse attendant-cub? = TRUE [report baseline-cost + (10.9 / 48) ] [ report baseline-cost + (2.6 / 48) ] ]  ;;;In MJ per tick
end

to-report  energy-cost-active
  let baseline-cost ( (0.010044 / 2) * body-mass)
  ifelse attendant-cub? = FALSE and attendant-yearling? = FALSE
[report baseline-cost]
[ifelse attendant-cub? = TRUE [report baseline-cost + (10.9 / 48) ] [ report baseline-cost + (2.6 / 48) ] ]  ;;;In MJ per tick
end

to-report nests-I-can-smell
    let x (patches in-radius 16)
  ask patch-here [set x other x]
     report (x with [count nests-here with [hidden? = FALSE] > ( (distance myself + 0.5) ^ 2) ] )
end

to-report random-sex
  ifelse random-float 1 > 0.5 [report "male"] [report "female"]
end

to-report random-binomial [n p]
     report sum n-values n [ifelse-value (p > random-float 1) [1] [0]]
end

to-report inverse-logit [ x ]
  report ( exp (x) / (1 + exp(x)) )
end
@#$#@#$#@
GRAPHICS-WINDOW
1118
61
4423
670
1000
175
1.647
1
14
1
1
1
0
0
0
1
-1000
1000
-175
175
1
1
1
ticks
30.0

BUTTON
1403
12
1467
45
go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

PLOT
1
14
310
223
Current Year Nest Status
Time
Count
0.0
100.0
0.0
100.0
true
true
"" "if day-of-year = 121 and hour-of-day = 0 [clear-plot]"
PENS
"active nests" 1.0 0 -13345367 true "" "plot active-nests"
"hatched nests" 1.0 0 -7500403 true "" "plot hatched-nests"
"failed nests" 1.0 0 -2674135 true "" "plot failed-nests"
"depredated nests " 1.0 0 -11085214 true "" "plot depredated-nests"

BUTTON
1122
10
1258
43
SetUp Landscape
setup-landscape
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1262
10
1395
43
Set Up Agents
setup-agents
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
1292
94
1384
139
NIL
hour-of-day
17
1
11

MONITOR
1197
94
1285
139
NIL
day-of-year
17
1
11

MONITOR
3
230
127
275
NIL
ice-breakup-date
0
1
11

MONITOR
5
330
127
375
mean nesting
mean-nesting-date
17
1
11

MONITOR
1136
94
1193
139
Year
year
17
1
11

PLOT
658
10
1013
189
Nest predation by year
Year
Number of nests
1.0
50.0
0.0
5000.0
true
false
"" ""
PENS
"depredated-nests" 1.0 0 -13840069 true "" "if year >= 2 [\nplotxy (item (year - 2) year-list) (item (year - 2) depredated-nests-by-year )\n]"

PLOT
283
237
645
364
Ice-breakup date by year
Year
Julian Date
1.0
50.0
130.0
200.0
false
false
"" ""
PENS
"Ice-breakup-date" 1.0 0 -16777216 true "" "if year >= 2 [ \nplotxy (item (year - 2) year-list) (item (year - 2) ice-breakup-by-year)]"

MONITOR
3
280
127
325
mean bear arrival
ice-breakup-date + 25
17
1
11

PLOT
287
364
645
488
freezeup date by year
Year
Julian Date
1.0
50.0
280.0
340.0
false
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "if year >= 2[\nplotxy (item (year - 2 ) year-list) (item (year - 2) ice-freezeup-by-year)\n]"

SWITCH
325
15
437
48
bears-on?
bears-on?
0
1
-1000

PLOT
650
192
1082
381
Eider popn size
Year
Count
1.0
50.0
0.0
10000.0
false
true
"" ""
PENS
"Adults" 1.0 0 -2674135 true "" "if year >= 2 [plotxy (item (year - 2) year-list) (item (year - 2) adults-by-year)]"
"Fledglings" 1.0 0 -817084 true "" "if year >= 2 [plotxy (item (year - 2) year-list) (item (year - 2) fledglings-by-year)]"
"Juveniles" 1.0 0 -13840069 true "" "if year >= 2 [plotxy (item (year - 2) year-list) (item (year - 2) juveniles-by-year)]"

PLOT
667
383
1120
533
Eider vital rates
Year
Percentage
1.0
50.0
0.0
100.0
false
true
"" ""
PENS
"Breeding pro" 1.0 0 -13345367 true "" "if year >= 2 [plotxy (item (year - 2) year-list) (item (year - 2) breeding-propensity-by-year * 100)]"
"Hatchling surv" 1.0 0 -7500403 true "" "if year >= 2 [plotxy (item (year - 2) year-list) (item (year - 2) hatchling-survival-by-year * 100)]"
"Clutch size (x 10)" 1.0 0 -2674135 true "" "if year >= 2 [plotxy (item (year - 2) year-list) (item (year - 2) clutch-size-by-year * 10)]"

SWITCH
325
53
448
86
bear-decline?
bear-decline?
1
1
-1000

SWITCH
323
155
594
188
clutch-size-climate-link?
clutch-size-climate-link?
0
1
-1000

SWITCH
323
118
594
151
breeding-propensity-climate-link?
breeding-propensity-climate-link?
0
1
-1000

@#$#@#$#@
This is a model of how polar bears and common eider ducks might respond to declines in Arctic sea ice.

For model details see the full ODD model description associated with this model at (XXX)


To run this model:

-make sure you have the rngs extension downloaded (https://github.com/cstaelin/RNGS-Extension)

-use Netlogo 5.3.1

Then,

1. Set desired model conditions
	-bears on?
	-bear decline?
	-clutch-size-climate-link?
	-breeding-propensity-climate-link?
2. Set Up Landscape (may take a bit)
3. Set Up Agents
4. Go
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

bear
false
15
Circle -1 true true 203 65 88
Circle -1 true true 137 106 120
Polygon -7500403 true false -82 120 -60 165 -45 165 -22 120
Circle -16777216 true false 248 111 16
Rectangle -1 true true 153 218 177 291
Polygon -1 true true 45 285 30 285 26 240 15 195 45 210
Rectangle -1 true true 60 221 80 292
Polygon -1 true true 194 285 209 285 218 243 239 210 194 210
Polygon -1 true true 279 82 288 102 305 96 297 80
Polygon -1 true true 220 78 211 98 194 92 202 76
Polygon -1 true true 240 135 255 135 105 120
Polygon -16777216 true false 254 129 265 129 254 137 243 128
Circle -16777216 true false 260 90 10
Circle -16777216 true false 235 90 10
Polygon -16777216 true false 287 91 292 88 294 94
Polygon -16777216 true false 208 83 206 91 216 90
Polygon -1 true true 91 85 143 92 161 100 184 108 205 115 214 208 204 217 180 225 154 233 132 232 109 229
Circle -1 true true 6 79 160

bird
false
0
Polygon -7500403 true true 135 165 90 270 120 300 180 300 210 270 165 165
Rectangle -7500403 true true 120 105 180 237
Polygon -7500403 true true 135 105 120 75 105 45 121 6 167 8 207 25 257 46 180 75 165 105
Circle -16777216 true false 128 21 42
Polygon -7500403 true true 163 116 194 92 212 86 230 86 250 90 265 98 279 111 290 126 296 143 298 158 298 166 296 183 286 204 272 219 259 227 235 240 241 223 250 207 251 192 245 180 232 168 216 162 200 162 186 166 175 173 171 180
Polygon -7500403 true true 137 116 106 92 88 86 70 86 50 90 35 98 21 111 10 126 4 143 2 158 2 166 4 183 14 204 28 219 41 227 65 240 59 223 50 207 49 192 55 180 68 168 84 162 100 162 114 166 125 173 129 180

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.3.1
@#$#@#$#@
setup
set grass? true
repeat 75 [ go ]
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="RunB" repetitions="100" runMetricsEveryStep="false">
    <setup>setup-landscape
setup-agents</setup>
    <go>go</go>
    <exitCondition>year &gt; 50</exitCondition>
    <metric>depredated-nests-by-year</metric>
    <metric>ice-breakup-by-year</metric>
    <metric>ice-freezeup-by-year</metric>
    <metric>mean-nest-ycor</metric>
    <metric>var-nests-per-cell</metric>
    <metric>nests-in-large-colonies</metric>
    <metric>incidence-abundance</metric>
    <metric>breeding-propensity-by-year</metric>
    <metric>clutch-size-by-year</metric>
    <metric>hatchling-survival-by-year</metric>
    <metric>nests-per-island-by-year</metric>
    <metric>adults-by-year</metric>
    <metric>fledglings-by-year</metric>
    <metric>juveniles-by-year</metric>
    <metric>mean-age-adults-by-year</metric>
    <metric>mean-age-adults-and-juveniles-by-year</metric>
    <enumeratedValueSet variable="bears-on?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CL-BP?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bear-decline?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="clutch-size-climate-link?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
