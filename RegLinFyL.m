%Script de prueba. NO TESTEADO AUN. 18/12/2020
%ENTRADA
%Leer archivo de entrada: lista de directorios \t PM(CH3X) \t Masa(C*)
%%Estructura archivo *.out para CH3X solo:
[file,PMads,Msol] = textread ('FILES.out', '%s %f %f'); %% PMads: peso molecular del CH3X; Msol: masa del Sólido C*
for i=1:length(file);
fileISOT = char(file(i));
partes = strsplit(fileISOT,{'Graphic_inputs/','_','.out'}) % Guarda las partes en un cell array
CH3X = partes{2}
Type = partes{3}
Ads = partes{4}
Temp = partes{5}
namePART = strsplit(fileISOT,{'Graphic_inputs/','.out'})
nameISOT = namePART{2}
isoterma = load(fileISOT);
Nads0 = isoterma(:,2);
item = min(find(Nads0>=10^-2));
Pressure = isoterma(item:end,1);
Nads = isoterma(item:end,2);
%FILTRO Nads>=10^-2 y CAMBIO UNIDADES (mmol/gC*)
    for j=1:length(Pressure);
        if Nads >= 10^(-2)
        nISOT = Nads / (6.023 * Msol(i)) 
        pISOT = Pressure 
        endif
    endfor
%LINEARIZACIÓN LANGMUIR
puntosL = length(nISOT)
while true
  %%%%%% Lectura de Nads desde el inicio hasta "puntosL"
  pISOTL = pISOT(1:puntosL,1);
  nISOTL = nISOT(1:puntosL,1);
  %%%%%%
yL = pISOTL ./ nISOTL
xL = pISOTL
Lang = polyfit(xL,yL,1)
FLang = Lang(2) + Lang(1) .* pISOTL
R2Lang = corr(yL,FLang)^2
fprintf('L0=%f L1=%f R2L=%f\n',Lang(2),Lang(1),R2Lang)
%Aceptación de linearización: si R^2>=0.95, si no rechazar y eliminar el último punto, repetir el cálculo
puntosL = puntosL - 1
    if R2Lang >= 0.9
    break
    end
NmaxL = 1 / Lang(1)
bL = Lang(1) / Lang(2)
%Guardo Valores (Lang(2),Lang(1),R2L,NmaxL,bL,puntosL) en archivo.
nameL = strcat(nameISOT,"L.txt")
delimiter = '\t'
Contenido = strjoin({nameL,nameISOT,num2str(Lang(2)),num2str(Lang(1)),num2str(R2Lang),num2str(NmaxL),num2str(bL),num2str(puntosL)},'\t');
dlmwrite(nameL,Contenido,delimiter,"-append");
% fprintf(fid, format, value, ...)
end
%
%LINEARIZACIÓN FREUNDLICH
puntosF = length(nISOT)
while true
  %%%%%% Lectura de Nads desde el inicio hasta "puntosF"
  pISOTF = pISOT(1:puntosF,1);
  nISOTF = nISOT(1:puntosF,1);
  %%%%%%
yF = log(nISOTF)
xF = log(pISOTF)
Freund = polyfit(xF,yF,1)
FFreund = Freund(2) + Freund(1) .* pISOTF
R2Freund = corr(yF,FFreund)^2
fprintf('F0=%f F1=%f R2F=%f\n',Freund(2),Freund(1),R2Freund)
%Aceptación de linearización: si R^2>=0.95, si no rechazar y eliminar el último punto, repetir el cálculo
puntosF = puntosF - 1
    if R2Freund >= 0.9
    break
    end
nF = 1 / Freund(1)
KF = 10^Freund(2)
%GUARDO VALORES (Freund(2),Freund(1),R2F,nF,KF,puntosF) en archivo.
nameF = strcat(nameISOT,"F.txt")
delimeter = '\t'
dlmwrite(nameF,nameISOT,Freund(2),Freund(1),R2Freund,nF,KF,puntosF,delimiter,"-append")
end
%PLOTEO: Nads [mmol/g(C*)] vs p [Pa]
yLm = NmaxL .* (pISOTL * bL ./ (1 .+ pISOTL * bL) )
yFm = KF .* ( pISOTF .^ (1 / nF) )
plot(xISOT,yISOT,'*','markersize',8,'color','black','linestyle','--','linewidth',1); hold on;
plot(xISOT,yLm,'^','markersize',6,'color',[0.98,0.65,0.1],'linestyle','-','linewidth',1); hold on;
plot(xISOT,yFm,'o','markersize',6,'color',[0.45,0.75,0.27],'linestyle','-','linewidth',1);
%title();
legend('Simulated','Langmuir','Freundlich','location','northwest');
xlabel('P [Pa]');
ytext = strcat(CH3X,' [mmol/g\_{C*}]')
ylabel(ytext);
end
%
