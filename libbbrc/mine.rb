require 'bbrc'

if ARGV.size < 2
  exit
end

smi = ARGV[0]
cla = ARGV[1]

MyFminer = Bbrc::Bbrc.new()

file = File.new(smi, "r")
while (line = file.gets)
  id = line.split[0].to_i
  smi = line.split[1].to_s
#  puts "'#{id}|#{smi}'"
  MyFminer.AddCompound(smi, id)
end
file.close

file = File.new(cla, "r")
while (line = file.gets)
  id = line.split[0].to_i
  act = line.split[line.split.size-1].to_f
  MyFminer.AddActivity(act, id)
end
file.close

# gather results for every root node in vector instead of immediate output
MyFminer.SetConsoleOut(false)
(0 .. MyFminer.GetNoRootNodes()-1).each do |j|
  result = MyFminer.MineRoot(j)
  result.each do |r|
    puts r
  end
end

