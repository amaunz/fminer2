require 'fminer'
 MyFminer = Fminer::Fminer.new()
 MyFminer.AddCompound("COC1=CC=C(C=C1)C2=NC(=C([NH]2)C3=CC=CC=C3)C4=CC=CC=C4" , 1)
 MyFminer.AddCompound("O=C1NC(=S)NC(=O)C1C(=O)NC2=CC=CC=C2" , 2)
    # ... continue adding compounds
 MyFminer.AddActivity(true, 1)
 MyFminer.AddActivity(false, 2)
    # ... continue adding activities (true for active, false for inactive)
 print MyFminer.GetNoCompounds()  
 puts " compounds"
 # Toy example: special settings for mining all fragments
 # use no significance constraint
 MyFminer.SetChisqSig(0) 
 # refine structures with support 1
 MyFminer.SetRefineSingles(true) 
 # gather results for every root node in vector instead of immediate output
 MyFminer.SetConsoleOut(false)
 (0 .. MyFminer.GetNoRootNodes()-1).each do |j|
    result = MyFminer.MineRoot(j)
    puts "Results"
    result.each do |res|
        puts res
   end
 end

