classdef PcacClass < handle
    %PCACCLASS wrapper for PCAC related classes/functions  
    
    properties
        opts %options struct containing parameters
        idMethod %string stating the ID method to be used
        mpcMethod %string stating the mpc/optimization method to be used
        idObj %object that computes theta
        optObj %object that computes control 
    end
    
    methods
        function obj = PcacClass(idMethod,mpcMethod,opts_, r_, specialFlag_ )
            %PCACCLASS Construct an instance of this class
            %   initialize parameters and the objects
            %r_ is the first commanded step
            obj.opts = opts_;
            
            %save the method to be used
            obj.idMethod = idMethod;
            obj.mpcMethod = mpcMethod;
            
            %set the id and opt obj
            %ID
            if strcmp(obj.idMethod, 'Linear_RLS')
                obj.idObj = PCAC.LinearRlsId(opts_,specialFlag_); %i.e. use LinearId from the PCAC package
            else
                error('Invalid ID method chosen')
            end
            %OPT
            if strcmp(obj.mpcMethod, 'QP')
                obj.optObj = PCAC.OptClass(opts_, r_,specialFlag_ );
            elseif strcmp(obj.mpcMethod, 'MATLAB')
                %Use matlab's mpc toolbox: needs a model to be initialized.
                %The obj will be initialized in the onestep method
            else
                error('Invalid MPC method chosen')
            end 
        end
        
        function [uOut,UOut] = oneStep(obj, u_, y_, r_ ,specialFlag_, UC)
            %oneStep computes one step of the PCAC update
            
            %update/get ID
            if strcmp(obj.idMethod, 'Linear_RLS')
                obj.idObj.updateCoeff(u_, y_)
                [ Fbar , Gbar ] = obj.idObj.getFG();
                [ Yout , Uout ] = obj.idObj.getIObuffers();
            end

            
            %update control
            if strcmp(obj.mpcMethod, 'QP')
                if specialFlag_
                    obj.optObj.updateUC(UC);
                end
                obj.optObj.updateControl( Fbar , Gbar , Yout , Uout , r_ ,specialFlag_)
                [uOut,UOut] = obj.optObj.getControl;
            end
            
        end
    end
end

