class RaxmlResultsParser

  attr_reader :data

  def initialize(file)
    @data = []
    parse(file)
  end

  def parse(file)
    f = File.open(file,'r'){|line|
      while line.gets
        @data << $_
        puts $_
      end
    }    
  end

end
