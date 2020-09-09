function Ai = SDHmatrix(a,alpha,q,d)
%SDHMATRIX Summary of this function goes here
%   Transfrom matrix in SDH coordinate
%   Output: 
%        Ai  --- transform matrix between 2 coordinate
%   Input:
%        a   --- lenghth along links 
%      alpha --- rotation angle 
%        q   --- anlge between links
%        d   --- length along z-axis

Ai = [cos(q), -sin(q)*cos(alpha), sin(q)*sin(alpha), a*cos(q);
      sin(q),  cos(q)*cos(alpha),-cos(q)*sin(alpha), a*sin(q);
      0     ,   sin(alpha)      ,  cos(alpha)       ,   d    ;
      0     ,       0           ,       0           ,   1     ];
  
end

