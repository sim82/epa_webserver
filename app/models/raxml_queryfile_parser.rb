class RaxmlQueryfileParser

  attr_reader :data, :valid_format, :error  

  def initialize(stream)
    @filename = stream.original_filename
    @data = stream.readlines
    @error = ""
    @valid_format = false
    check_format
  end

  def check_format
    i = 0
    while i < @data.size
      if @data[i] =~ /^>/ #fasta format
        j = i+1
        while @data[j]=~/^[a-zA-Z]+$/ && j < @data.size 
          j = j+1
        end
      else
        @error = "Queryfile is not in FASTA format!! \n ParserError :: #{@filename} line: #{i} => #{@data[i]}"
        return
      end
      i = j
    end
    @valid_format = true
  end
end
