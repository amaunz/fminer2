require './bbrc'
 MyFminer = Bbrc::Bbrc.new()
 MyFminer.SetConsoleOut(false)

 # ---- SETTINGS. ----
 #
 # Toy example: special settings for mining all fragments
 # see http://opentox.github.com/algorithm/2012/05/02/use-case-table-for-fminer/
 MyFminer.SetChisqSig(0) 
 MyFminer.SetBackbone(0)
 MyFminer.SetDynamicUpperBound(0)

 # Refine structures with support according to weights!
 # Try values between 2.0 and 8.0
 MyFminer.SetMinfreq(3.1) # IMPORTANT! TRY DIFFERENT VALUES!
 # Try 0.9, 1.0, 1.1
 # Try 2.9, 3.0, 3.1
 # Try 6.9, 7.0, 7.1
 # for weights of 4.0, 2.0, 1.0, as below

 # ---- ADD DATA BELOW. DO NOT USE "Set*" METHODS BELOW THIS LINE ----

 # Add compounds below.
 MyFminer.AddCompound("C-C" , 1)
 MyFminer.AddCompound("C-C-C" , 2)
 MyFminer.AddCompound("C-C-C-C" , 3)
 # ... continue adding compounds
 MyFminer.AddActivity(0.0, 1)
 MyFminer.AddActivity(0.0, 2)
 MyFminer.AddActivity(0.0, 3)
 # ... continue adding activities
 # ... treated as categories unless MyFminer.SetRegression(1) is used
 MyFminer.AddWeight(4.0, 1);
 MyFminer.AddWeight(2.0, 2);
 MyFminer.AddWeight(1.0, 3);
 # ... continue adding weights (positive floats >= 1.0)
 # ... weights are rounded down, e.g. 1.2 => 1.0 and 1.8 => 1.0
 
 # ---- START MINING. ----

 # gather results for every root node in vector instead of immediate output
 (0 .. MyFminer.GetNoRootNodes()-1).each do |j|
    result = MyFminer.MineRoot(j)
    result.each do |res|
        puts res
   end
 end

