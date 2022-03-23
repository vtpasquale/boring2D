classdef Matrix3D < double
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    methods
        function obj = Matrix3D(valueIn)
        	if ~isa(valueIn,'double'); error('Incorrect input type'); end
            obj = obj@double(valueIn);
        end
        function C = mtimes(A,B)
            % A : (a x c x Z)
            % B : (c x b x Z)
            Ap = permute(A,[2,1,4,3]); % (c x a x 1 x Z)
            Bp = permute(B,[1,4,2,3]); % (c x 1 x b x Z)
            C = Ap .* Bp;              % (c x a x b x Z)
            C = sum(C,1);              % (1 x a x b x Z)
            C = permute(C,[2,3,4,1]);  % (a x b x Z)
        end
        function B = transpose(A)
            B = permute(A,[2,1,3]);
        end
%         function C = benchTimes(A,B)
%             
%             tic
%             C = mtimes(A,B);
%             tocPermute = toc
%             
%             tic 
%             C1 = mTimesBigLoop(A,B);
%             tocBigLoop = toc
%             
%             tic 
%             C2 = mTimesDoubleLoop(A,B);
%             tocDoubleLoop = toc
%             
%             max(abs(C2(:)-C(:)))
%             max(abs(C1(:)-C(:)))
%         end
%         function C = mTimesBigLoop(A,B)
%             % A : (a x c x Z)
%             % B : (c x b x Z)
%             A = double(A);
%             B = double(B);
%             nA = size(A);
%             nB = size(B);            
%             if nA(3) ~= nB(3); error('nA(3) ~= nB(3)'); end
%             if nA(2) ~= nB(1); error('Inner matrix dimensions must agree'); end
%             C = zeros(nA(1),nB(2),nA(3));
%             for k = 1:nA(3)
%                 C(:,:,k) = A(:,:,k)*B(:,:,k);
%             end
%         end
%         function C = mTimesDoubleLoop(A,B)
%             % A : (a x c x Z)
%             % B : (c x b x Z)
%             nA = size(A);
%             nB = size(B);            
%             if nA(3) ~= nB(3); error('nA(3) ~= nB(3)'); end
%             if nA(2) ~= nB(1); error('Inner matrix dimensions must agree'); end
%             C = zeros(nA(1),nB(2),nA(3));
%             for row = 1:nA(1)
%                 for col = 1:nB(2)
%                     C(row,col,:) = sum( A(row,:,:).*(B(:,col) ).',2);
%                 end
%             end
%         end
    end
    
end

