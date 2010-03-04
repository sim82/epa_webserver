module ActiveRecord
  module Validators
    module ClassMethods
      
      def validates_heuristic(*attr_names)
        configuration = {
          :message => "Invalid value entered",
          :allow_blank => true,
          :on => :save}

        configuration.update(attr_names.pop) if attr_names.last.is_a?(Hash)
        
      end

    end
  end
end
