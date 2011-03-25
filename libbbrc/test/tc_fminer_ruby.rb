require "test/unit"
require 'digest/md5'
require 'yaml'

begin
  require File.dirname(__FILE__) +  "/fminer.rb"
rescue Exception
  puts File.new(__FILE__).path + ": fminer.rb not found!"
  exit false
end


class TestSimpleNumber < Test::Unit::TestCase

    
  def test_ruby_fminer
    @smi_file=File.dirname(__FILE__) + "/hamster_carcinogenicity.smi"
    @class_file=File.dirname(__FILE__) + "/hamster_carcinogenicity.class"
    @md5_yaml_file=File.dirname(__FILE__) + "/fminer_default_md5.yaml"

    #assert true
    #puts @smi_file
    output=run_fminer(@smi_file, @class_file, 2)
    actual_md5=Digest::MD5.hexdigest(output)
    # puts actual_md5
    config = YAML.load_file(@md5_yaml_file)
    expected_md5=config['default']
    # puts expected_md5
    assert_equal(actual_md5, expected_md5)

  end

end
