function [C1, m] = tubebankpar(Sld, Std)


%parameters to determine the nusselt # from aligned tube bank
%From Incropera, pg. 438

C1par = [ .348 .275 .1 .0633;
          .367 .250 .101 .0678;
          .418 .299 .229 .198;
          .290 .357 .374 .286];
      
mpar =  [ .592 .608 .704 .752;
          .586 .620 .702 .744;
          .570 .602 .632 .648;
          .601 .548 .581 .608];
if Sld == 1.25 && Std == 1.25
    C1 = C1par(1,1);
    m  = mpar(1,1);
elseif Sld == 1.5 && Std == 1.25
    C1 = C1par(2,1);
    m  = mpar(2,1);
elseif Sld == 2 && Std == 1.25
    C1 = C1par(3,1);
    m  = mpar(3,1); 
elseif Sld == 3 && Std == 1.25
    C1 = C1par(4,1);
    m  = mpar(4,1);
    
elseif Sld == 1.25 && Std == 1.5
    C1 = C1par(1,2);
    m  = mpar(1,2);
elseif Sld == 1.25 && Std == 2
    C1 = C1par(1,3);
    m  = mpar(1,3);
elseif Sld == 1.25 && Std == 3
    C1 = C1par(1,4);
    m  = mpar(1,4);

elseif Sld == 1.5 && Std == 1.5
    C1 = C1par(2,2);
    m  = mpar(2,2);
elseif Sld == 1.5 && Std == 2
    C1 = C1par(2,3);
    m  = mpar(2,3);
elseif Sld == 1.5 && Std == 3
    C1 = C1par(2,4);
    m  = mpar(2,4);

elseif Sld == 2 && Std == 1.5
    C1 = C1par(3,2);
    m  = mpar(3,2);
elseif Sld == 2 && Std == 2
    C1 = C1par(3,3);
    m  = mpar(3,3);
elseif Sld == 2 && Std == 3
    C1 = C1par(3,4);
    m  = mpar(3,4);
      
elseif Sld == 3 && Std == 1.5
    C1 = C1par(4,2);
    m  = mpar(4,2);
elseif Sld == 3 && Std == 2
    C1 = C1par(4,3);
    m  = mpar(4,3);
elseif Sld == 3 && Std == 3
    C1 = C1par(4,4);
    m  = mpar(4,4);
    
else
    disp('Tube bank is out of parameter specs, change Sld and/or Std')
    
    
end
      
 end