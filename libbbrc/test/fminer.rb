# re-implementation of fminer in ruby with integrated last-utils
# args: smi, class 
# TODO mode, minfreq

begin
  require File.dirname(__FILE__) +  "/../bbrc.so"
rescue
  puts "last not found!"
  exit false
end

def read_tab_file(f)
  fc = []
  File.open(f).each do |line|
    line.chomp!
    fc << line.split("\t")
  end
  fc
end


def run_fminer(smi_file, class_file, min_freq)
  myFminer = Bbrc::Bbrc.new()
  # Adjust settings
  myFminer.SetConsoleOut(false)
  myFminer.SetMinfreq(min_freq.to_i)

  # Add compounds below. IMPORTANT! DO NOT CHANGE SETTINGS AFTER ADDING COMPOUNDS!
  smi_h = {} 
  smi=read_tab_file(smi_file)
  smi.each do |s|
    myFminer.AddCompound(s[1].to_s , s[0].to_i)
    smi_h[s[0].to_i] = s[1].to_s
  end

  cla=read_tab_file(class_file)
  all_hash = Hash.new
  cla.each do |c|
    id=c[0].to_i
    activity=c[2].to_f
    myFminer.AddActivity(activity , id)
    all_hash[id]=activity
  end

  # gather results for every root node in vector instead of immediate output
  result_str = ""
  (0 .. myFminer.GetNoRootNodes()-1).each do |j|
    result = myFminer.MineRoot(j)
    result.each do |res|
      result_str+=res+"\n"
    end
  end
  result_str

end
