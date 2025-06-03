PRO TAD

  path_img = '...' ;one image with several bands of a spectral index; each one representing one instant   
  path_out =  '...' ;output image
  xi = 0.1              ;method's parameter
  ;--------------------------------------------------------------------------

  img = read_tiff(path_img)

  nb = n_elements(img[*,0,0])
  nc = n_elements(img[0,*,0])
  nl = n_elements(img[0,0,*])

  transfImage = fltarr(nc,nl)

  t0 = systime(/seconds)

  ;Data transformation
  W = fltarr(nb) + 1
  W *= 1/norm(W)
  WT = transpose(W)
  nW = norm(W)

  for i = 0, nc-1 do begin
    for j = 0, nl-1 do begin
      projEscalar = ( img[*,i,j] ## WT ) / nW
      proj =  W # projEscalar
      d = norm(img[*,i,j] - proj)

      transfImage[i,j] = d
    endfor
  endfor

  ;Thresholding according Negri et al (2022)
  posGTMu = where(transfImage ge mean(transfImage) - xi*stddev(transfImage) )
  __transfImageVector = transfImage[posGTMu]
  resBB = best_bouding_box(__transfImageVector)
  statThres_sup = resBB.lims[1]
  tau = statThres_sup
  posc = where(transfImage gt tau)
  imTau = transfImage[*,*] * 0
  imTau[posc] = 1

  t1 = systime(/seconds)

  ;Dimension adjust
  Path_GEO = path_img
  Result = QUERY_TIFF(Path_GEO, Info, GEOTIFF=geoVar)
  dimRef = info.DIMENSIONS
  tauImage_CG = CONGRID(imTau, dimRef[0], dimRef[1])

  ;Saving the output--------------------------------------------
  write_tiff, path_out, geotiff=geoVar, tauImage_CG
  ;--------------------------------------------------------------

  print, 'Procesing time: ', t1-t0
  print, 'End of process...'
END



;------------------------------
function best_bouding_box, data

  k = 1
  seed = 1234567L
  rep = 10000L

  BS = 0.5*IQR(data)*(N_ELEMENTS(data)^(-1.0/3.0)) ;freedman-diaconis rule
  y = HISTOGRAM(data, BINSIZE=BS, LOCATIONS = x)

  ;Remove extreme values---------------
  z = total(y,/cumulative)/total(y)
  posInf = where(z lt 0.01) & if posInf[0] eq -1 then posInf[0] = 0
  posSup = where(z gt 0.99) & if posSup[0] eq -1 then posSup[0] = n_elements(z)-1
  y = y[posInf[n_elements(posinf)-1]:posSup[0]]
  x = x[posInf[n_elements(posinf)-1]:posSup[0]]
  ;------------------------------------

  nx = n_elements(x)
  dx = x[1]-x[0]

  ;sentinelPos = 10D^100
  sentinelPos = 0
  sentinelCost = (max(y) * (x[n_elements(x)-1] - x[0])) - total(y*dx)  ;sem divisão...
  ;vecStructBox = [ptr_new({lims: [min(x), max(x)], minCost: sentinelCost, bestConf: [0,n_elements(x)-1]})]

  cost = dblarr(rep)
  conf = lonarr(k+2,rep)

  for r = 0, rep-1 do begin

    rnd = sort(randomu(seed,nx))
    conf[*,r] = [0, rnd(sort(rnd[0:k-1])), nx-1]  ;represents the partition configuration

    penaltyDiag = 0
    for i = 0, k do begin
      BOX = max(y[conf[i,r]:conf[i+1,r]]) * (x[conf[i+1,r]] - x[conf[i,r]])
      AUF = total( y[conf[i,r]:conf[i+1,r]]*dx )
      cost[r] += (BOX - AUF)
    endfor

  endfor

  minPos = where(cost eq min(cost))
  bestConf = conf[*,minPos[0]]
  lims = x[bestConf]

  midPoints = fltarr(n_elements(bestConf)-1)
  for i = 0, n_elements(bestConf)-2 do midPoints[i] = (bestConf[i] + bestConf[i+1])*0.5

  return, {bestConf: bestConf, lims: lims, midPoints: midPoints} 
end


;------------------------
FUNCTION IQR, Img

  sortData = Img[SORT(Img)]
  ind = N_ELEMENTS(Img)/2

  IF N_ELEMENTS(sortData) MOD 2 EQ 0 THEN BEGIN
    lower = sortData[0:ind-1]
    higher = sortData[ind:N_ELEMENTS(Img)-1]
  ENDIF ELSE BEGIN
    lower = sortData[0:ind]
    higher = sortData[ind:N_ELEMENTS(Img)-1]
  ENDELSE

  q25 = MEDIAN(lower, /EVEN)
  q75 = MEDIAN(higher, /EVEN)

  Return, q75-q25
END