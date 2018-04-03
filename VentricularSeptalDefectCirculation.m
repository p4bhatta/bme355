classdef VentricularSeptalDefectCirculation < handle
    % A Non-linear State-Space Model of Model of systemic and pulmonary
    % circulation of a heart with ventricular septal defect
    
    properties (SetAccess = private)
        HR; %heart rate
        tc; %contraction time (seconds per beat)
        Tmax; 
    end
    
    properties (Access = public)
        EmaxL;
        EmaxR;
        Emin;
        
        nonSlackBloodVolume = 210; % (ml)
        
        %Resistances for the 
        R1 = .5; % between .5 and 2
        R2 = .005;
        R3 = .001;
        
        R4 = .5;
        R5 = .005; % between .5 and 2
        R6 = .001;
        
        R7 = .0398;
        R8 = .0398;

        C1= 4.4;
        C3 = 5.55;
        C5 = 1.33;
        
        L1 = .0005;
        L2 = .0005;
    end
    
    methods (Access = public)
        function C = VentricularSeptalDefectCirculation(HR, EmaxL, EmaxR, Emin) 
            % HR: heart rate (beats per minute)
            % Emax: maximum elastance
            % Emin: minimum elastance
            
            C.setHeartRate(HR);
            C.EmaxL = EmaxL;
            C.EmaxR = EmaxR;
            C.Emin = Emin;
        end
        
        function setHeartRate(C, HR)
            % beats per minute
            C.HR = HR;
            C.tc = 60/HR; %seconds per beat
            C.Tmax = .2+.15*C.tc; % contraction time
        end
        
        function dx = getDerivative(C, t, x, heartState)
            
            % x1 - right atrial pressure
            % x2 - right ventricular pressure
            % x3 - left atrial pressure
            % x4 - left ventricular pressure
            % x5 - arterial pressure
            % x6 - pulmonary artery flow
            % x7 - aortic flow

            % t: time
            % x: state variables [ventricular pressure; atrial pressure; arterial pressure; aortic flow]
            % dx: time derivatives of state variables
            % heartState: 0 - healthy heart, 1 - moderate septal defect, 2 - severe septal defect
            
            pulmonaryFlow = x(6);
            aorticFlow = x(7);
            t
            if x(1) > x(2) || x(3) > x(4) %flow == 0 
                disp("1 - Filling");
                A = filling(C, t, heartState);
            elseif (x(4) > x(5)) || (pulmonaryFlow > 0) || (aorticFlow > 0)%FINISHED
                disp("2 - Ejection");
                A = ejection(C, t, heartState);
            else %flow == 0
                disp("0 - Isovolumic");
                A = isovolumic(C, t, heartState);
            end
            dx = A*(x);
        end
        
        function A = isovolumic(C, t, heartState)
            % Produces the isovolumic A matrix for moderate Ventricular
            % Septal Defect (VSD)
            % t: time
            % isHealthy: 0 or 1 specifying if it is moderate or severe VSD

            %elastance of the left side of the heart
            [elR,elL] = elastance(C,t);
            delR_dt = elastanceFiniteDifference(C, t);
            
            %elastance of the right side of the heart
            delL_dt = elastanceFiniteDifference(C, t);
            
            C2 = 1/elL;
            C4 = 1/elR;
            
            if (heartState == 0) 
                A = [1/(C.C1*C.R8) 0 0 0 -1/(C.C1*C.R8) 0 0; %x1'
                0 0 0 0 0 0 0; %x2'
                0 0 0 0 0 0 0; %x3'
                0 0 0 0 0 0 0; %x4'
                1/(C.C5*C.R8) 0 0 0 -1/(C.C5*C.R8) 0 0; %x5'
                0 0 0 0 0 0 0; %x6'
                0 0 0 0 0 0 0]; %x7
            elseif (heartState == 1)
                A = [-1/(C.C1*C.R8) 0 0 0 1/(C.C1*C.R8) 0 0; %x1'
                0 -1/(C2*C.R7) 0 1/(C2*C.R7) 0 0 0; %x2'
                0 0 0 0 0 0 0; %x3'
                0 -1/(C4*C.R7) 0 1/(C4*C.R7) 0 0 0; %x4'
                -1/(C.C5*C.R8) 0 0 0 1/(C.C5*C.R8) 0 0; %x5'
                0 0 0 0 0 0 0; %x6'
                0 0 0 0 0 0 0]; %x7'
            else
                A = [1/(C.C1*C.R8) 0 0 0 -1/(C.C1*C.R8) 0 0; %x1'
                0 (elR/C.R7-delR_dt/elR) 0 -elR/C.R7 0 0 0; %x2'
                0 0 0 0 0 0 0; %x3'
                0 elL/C.R7 0 (-elL/C.R7-delL_dt/elL) 0 0 0; %x4'
                1/(C.C5*C.R8) 0 0 0 -1/(C.C5*C.R8) 0 0; %x5'
                0 0 0 0 0 0 0; %x6'
                0 0 0 0 0 0 0]; %x7'
            end
        end
        
        function A = filling(C, t, heartState)
            % Produces the isovolumic A matrix for moderate Ventricular
            % Septal Defect (VSD)
            % t: time
            % isSevere: 0 or 1 specifying if it is moderate or severe VSD

            % set resistor 7 
            if (heartState == 2) 
                res7 = 0; 
            else
                res7 = C.R7;
            end
            %elastance of the left side of the heart
            [elR,elL] = elastance(C,t);
            delR_dt = elastanceFiniteDifference(C, t);

            %elastance of the right side of the heart
            delL_dt = elastanceFiniteDifference(C, t);
            
            if (heartState == 0)
                A = [(1/C.R1 + 1/C.R2)/C.C1 -1/(C.C1*C.R8) 0 0 -1/(C.C1*C.R8) 0 0; %x1'
                    elR/C.R1 (delR_dt/elR - elR/C.R1) 0 0 0 0 0; %x2'
                    0 0 1/(C.C3*C.R4) -1/(C.C3*C.R4) 0 0 0; %x3'
                    0 0 -elL/C.R4 (elL/C.R4+delL_dt/elL) 0 0 0; %x4'
                    1/(C.C5*C.R8) 0 0 0 -1/(C.C5*C.R8) 0 0; %x5'
                    0 0 0 0 0 0 0; %x6'
                    0 0 0 0 0 0 0]; %x7'
            else
                A = [(1/(C.R8*C.C1)+1/(C.R1*C.C1)) -1/(C.C1*C.R1) 0 0 -1/(C.R8*C.C1) 0 0; %x1'
                    elR/C.R1 (delR_dt/elR+elR/C.R7-elR/C.R1) 0 -elR/C.R7 0 0 0; %x2'
                    0 0 1/(C.C3*C.R4) -1/(C.C3*C.R4) 0 0 0; %x3'
                    0 elL/C.R7 -elL/(C.R4) (elL/C.R4-elL/C.R7+delL_dt/elL) 0 0 0; %x4'
                    1/(C.C5*C.R8) 0 0 0 -1/(C.C5*C.R8) 0 0; %x5'
                    0 0 0 0 0 0 0; %x6'
                    0 0 0 0 0 0 0]; %x7'
            end
        end
        
        function A = ejection(C, t, heartState)
            % Produces the ejection A matrix for Ventricular
            % Septal Defect
            % t: time
            % isSevere: 0 or 1 specifying if it is moderate or severe VSD

            % set resistor 7
            if (heartState == 2) 
                res7 = 0;
            else
                res7 = C.R7;
            end

            %elastance of the left side of the heart
            [elR,elL] = elastance(C,t);
            delR_dt = elastanceFiniteDifference(C, t);

            %elastance of the right side of the heart
            delL_dt = elastanceFiniteDifference(C, t);
            
%             A = [1/(C.C1*C.R8) 0 0 0 1/(C.C1*C.R8) 0 0; %x1'
%                  0 (-1/C.R7-delR_dt) 0 1/C.R7 0 1 0; %x2'
%                  0 0 0 0 0 -1/C.C3 0; %x3'
%                  0 1/C.R7 0 (-delR_dt-1/C.R7) 0 0 1; %x4'
%                  1/(C.C5*C.R8) 0 0 0 -1/(C.R8*C.C5) 0 -1/C.C5; %x5'
%                  0 -1/C.L1 1/C.L1 0 0 -(C.R2+C.R3)/C.L1 0; %x6'
%                  0 0 0 -1/C.L2 1/C.L2 0 -(C.R5+C.R6)/C.L2]; %x7'
%             
%             A = [1/(C.C1*C.R8) 0 0 0 -1/(C.C1*C.R8) 0 0; %x1'
%                  0 (delR_dt/elR - elR/C.R7) 0 elR/C.R7 0 -elR 0; %x2'
%                  0 0 0 0 0 1/C.C3 0; %x3'
%                  0 elL/C.R7 0 (elL/delL_dt^2 - elL/C.R7)  0 0 elL; %x4'
%                  -1/(C.C5*C.R8) 0 0 0 1/(C.R8*C.C5) 0 1/C.C5; %x5'
%                  0 -1/C.L1 -1/C.L1 0 0 -(C.R2+C.R3)/C.L1 0; %x6'
%                  0 0 0 -1/C.L2 -1/C.L2 0 -(C.R5+C.R6)/C.L2]; %x7'
            if (heartState == 0)
              A = [1/(C.C1*C.R8) 0 0 0 -1/(C.C1*C.R8) 0 0; %x1'
                 0 delR_dt/elR 0 -0 0 elR 0; %x2'
                 0 0 0 0 0 1/C.C3 0; %x3'
                 0 0 0 delL_dt/elL  0 0 elL; %x4'
                 -1/(C.C5*C.R8) 0 0 0 1/(C.R8*C.C5) 0 1/C.C5; %x5'
                 0 -1/C.L1 1/C.L1 0 0 -(C.R2+C.R3)/C.L1 0; %x6'
                 0 0 0 -1/C.L2 1/C.L2 0 -(C.R5+C.R6)/C.L2]; %x7'
            else
                A = [1/(C.C1*C.R8) 0 0 0 -1/(C.C1*C.R8) 0 0; %x1'
                 0 (delR_dt/elR+elR/C.R7) 0 -elR/C.R7 0 elR 0; %x2'
                 0 0 0 0 0 1/C.C3 0; %x3'
                 0 -elL/C.R7 0 (delL_dt/elL+elL/C.R7)  0 0 elL; %x4'
                 -1/(C.C5*C.R8) 0 0 0 1/(C.R8*C.C5) 0 1/C.C5; %x5'
                 0 -1/C.L1 1/C.L1 0 0 -(C.R2+C.R3)/C.L1 0; %x6'
                 0 0 0 -1/C.L2 1/C.L2 0 -(C.R5+C.R6)/C.L2]; %x7'
            end
        end
        
        function [elR,elL] = elastance(C, t)
            % t: time
            % result: time-varying elastance
            
            % C.Tmax is the length of time of the ventricular contraction 
            % tn is time normalized to C.Tmax
            tn = rem(t,C.tc)/C.Tmax;
            neg = find(tn < 0); %this line and the next are for generality 
            tn(neg) = tn(neg) + C.Tmax; % in case you give a t < 0
            En = 1.55 * (tn/.7).^1.9 ./ (1 + (tn/.7).^1.9) ./ (1 + (tn/1.17).^21.9); %En is normalized elastance (between 0 and 1)
            elL = (C.EmaxL-C.Emin)*En+C.Emin; % this rescales the normalized elastance to the range (Emin,Emax)
            elR = (C.EmaxR-C.Emin)*En+C.Emin; 
        end
        
        function result = elastanceFiniteDifference(C, t)
            % t: time
            % result: finite-difference approximation of time derivative of 
            % time-varying elastance            

            dt = .0001;
            forwardTime = t + dt;
            backwardTime = max(0, t - dt); % small negative times are wrapped to end of cycle
            forward = elastance(C, forwardTime);
            backward = elastance(C, backwardTime);
            result = (forward - backward) ./ (forwardTime - backwardTime);
        end
        
        function simulate(C, simulationTime, heartState)
            % simulationTime: # seconds to simulate
            % heartState: 0 - healthy heart, 1 - moderate septal defect, 2 - severe septal defect
            % We put all the blood pressure in the atria as an
            % initial condition. Note that we can't get the total blood
            % volume by multiplying by C2, because we're missing the
            % pulmonary loop. 
            
            initialState = [C.nonSlackBloodVolume/C.C1; 0; C.nonSlackBloodVolume/C.C3; 0; 0; 0; 0]; 
            ofun = @(t,x) C.getDerivative(t,x, heartState); % wrapper function because ode45 expects a function rather than a method
            [time, state] = ode45(ofun, [0 simulationTime], initialState);

            figure(1);
            hold on;
            plot(time, state(:, 1), 'r');
            plot(time, state(:, 2), 'g');
            plot(time, state(:, 3), 'b');
            title('Right Side Pressures');
            legend('atrial pressure', 'ventricular pressure', 'pulmonary pressure');
            set(gca, 'FontSize', 18);
            xlabel('Time (s)');
            ylabel('Pressure (mmHg)');
            
            figure(2);
            hold on;
            plot(time, state(:, 3), 'r');
            plot(time, state(:, 4), 'g');
            plot(time, state(:, 7) + state(:, 4)*C.R6, 'b');
            title('Left Side Pressures');
            legend('atrial pressure', 'ventricular pressure', 'aortic pressure');
            set(gca, 'FontSize', 18);
            xlabel('Time (s)');
            ylabel('Pressure (mmHg)');
        end
        
    end
        
end