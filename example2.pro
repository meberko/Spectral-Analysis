
; To run this program
;
; IDL> @example.pro
;
; JaeSub Hong, 2003-2005, version 1.5
; Please report any problem or suggestion at jaesub@head.cfa.harvard.edu
;

range = [0.3,8.0]                           ; 0.3 - 8.0 keV
frac  = [0.25, 0.33, 0.50, 0.67, 0.75]
nofixerror = 0                              ; set 0 to fix error and 1 to not fix error

readcol,'/Users/chandra/quantile_1.7/sources_lt_25_nc_gt_100.txt',srcnos_lt_25,format='A' ; get list of all source numbers we will be analyzing
print, srcnos_lt_25
n_src_lt_25=n_elements(srcnos_lt_25)
QDx_array_lt_25=fltarr(n_src_lt_25)
QDy_array_lt_25=fltarr(n_src_lt_25)
QDxerr_array_lt_25=fltarr(n_src_lt_25, 2)
QDyerr_array_lt_25=fltarr(n_src_lt_25, 2)

readcol,'/Users/chandra/quantile_1.7/sources_gt_25_nc_gt_100.txt',srcnos_gt_25,format='A' ; get list of all source numbers we will be analyzing
print, srcnos_gt_25
n_src_gt_25=n_elements(srcnos_gt_25)
QDx_array_gt_25=fltarr(n_src_gt_25)
QDy_array_gt_25=fltarr(n_src_gt_25)
QDxerr_array_gt_25=fltarr(n_src_gt_25, 2)
QDyerr_array_gt_25=fltarr(n_src_gt_25, 2)

j=0

for i=0, 33 do begin

;---------------------------------------------------------
; example of simple source photons
   readcol,'/Users/chandra/ae_run/extract/point_sources/quantile_inputs/'+srcnos_lt_25[i]+'_src.txt',src,format='f'
   Ex = quantile(frac, src, range=range, err_Ex=err_Ex, $
                 Qx=Qx, err_Qx=err_Qx, $
                 nofixerror=nofixerror, $
                 QDx=QDx, QDy=QDy, $
                 err_QDx=err_QDx, err_QDy=err_QDy)

;   print,'----------------------'
;   print,'src   ',n_elements(src)
;   print,'frac  ',frac
;   print,'Ex    ',Ex
;   print,'err_Ex',err_Ex
;   print,'Qx    ',Qx
;   print,'err_Qx',err_Qx
;   print,'QDx   ',QDx, err_QDx
;   print,'QDy   ',QDy, err_QDy

;---------------------------------------------------------
; example of bkgnd subtraction
   readcol,'/Users/chandra/ae_run/extract/point_sources/quantile_inputs/'+srcnos_lt_25[i]+'_bkg.txt', bkg, format='f'

;   ratio = 0.2 ; default set by JS
   ratio = 0.02
   Ex = quantile(frac, src, bkg, ratio=ratio, $
                 range=range, err_Ex=err_Ex, $
                 Qx=Qx, err_Qx=err_Qx, $
                 nofixerror=nofixerror, $
                 QDx=QDx, QDy=QDy, $
                 err_QDx=err_QDx, err_QDy=err_QDy)

   print,'----------------------'
   print,'src   ',n_elements(src)
   print,'bkg   ',n_elements(bkg)
   print,'net   ',n_elements(src)-ratio*n_elements(bkg)
   print,'ratio ',ratio
   print,'frac  ',frac
   print,'Ex    ',Ex
   print,'err_Ex',err_Ex
   print,'Qx    ',Qx
   print,'err_Qx',err_Qx
   print,'QDx   ',QDx, err_QDx
   print,'QDy   ',QDy, err_QDy

   QDx_array_lt_25(j)=QDx
   QDy_array_lt_25(j)=QDy
   QDxerr_array_lt_25(j, *)=err_QDx
   QDyerr_array_lt_25(j, *)=err_QDy
   j = j+1
endfor

;---------------------------------------------------------
; example of plotting QCCD
; to compile the plotting routines
qdplot


set_plot,'ps'
device,/portrait,filename='quantile.eps'


!P.charsize=1.0
!P.thick=1.5
!P.charthick=1.3

; rainbow color scheme
loadct,39

; set the plot window
xr    =[-0.25,0.3]
yr    =[1.4,2.2]
xtick =['0.6','1','2','3','4','5','  6 keV','7'] ; top x-axis
qd_setwin, xr=xr, yr=yr, range=range, xtick=xtick, stu=0.5

gamma = "103B
xyouts, -0.057, 1.42, '!4' + String(gamma) + '!X = 0', COLOR=250
xyouts, -0.132, 1.45, '!4' + String(gamma) + '!X = 1', COLOR=250
xyouts, -0.19, 1.51, '!4' + String(gamma) + '!X = 2', COLOR=250
xyouts, -0.22, 1.615, '!4' + String(gamma) + '!X = 3', COLOR=250
xyouts, -0.235, 1.717, '!4' + String(gamma) + '!X = 4', COLOR=250
xyouts, -0.23, 1.53, '4', COLOR=80, ORIENTATION=340

; read grid data
file_grid = '/Users/chandra/quantile_1.7/example_powerlaw_grid.rdb' ; use grid.pl to generate it
grid=qd_readgrid(file_grid, comment=comment)

; label grid
align=[1.0,0.5]	  ; label alignment
ori  =[-20,0]	  ; orientation
par  =[0.01,4]	  ; location
xo   =[0.0,0.0]	  ; x offset
yo   =[0.0,-0.05] ; y offset
col = [80,250]	  ; color setting
prefix1= ['','','','','','N!dH!n(x10!u22!n)=']
prefix2= '!4C!3='
qd_labelgrid, par, grid, comment=comment, $
	prefix1=prefix1, prefix2=prefix2, $
	align=align, ori=ori, col=col, xo=xo, yo=yo

; plot grid
qd_plotgrid, grid, comment=comment, col=col

j=0

for i=0, 57 do begin

;---------------------------------------------------------
; example of simple source photons
   readcol,'/Users/chandra/ae_run/extract/point_sources/quantile_inputs/'+srcnos_gt_25[i]+'_src.txt',src,format='f'
   Ex = quantile(frac, src, range=range, err_Ex=err_Ex, $
                 Qx=Qx, err_Qx=err_Qx, $
                 nofixerror=nofixerror, $
                 QDx=QDx, QDy=QDy, $
                 err_QDx=err_QDx, err_QDy=err_QDy)

;   print,'----------------------'
;   print,'src   ',n_elements(src)
;   print,'frac  ',frac
;   print,'Ex    ',Ex
;   print,'err_Ex',err_Ex
;   print,'Qx    ',Qx
;   print,'err_Qx',err_Qx
;   print,'QDx   ',QDx, err_QDx
;   print,'QDy   ',QDy, err_QDy

;---------------------------------------------------------
; example of bkgnd subtraction
   readcol,'/Users/chandra/ae_run/extract/point_sources/quantile_inputs/'+srcnos_gt_25[i]+'_bkg.txt', bkg, format='f'

;   ratio = 0.2 ; default set by JS
   ratio = 0.02
   Ex = quantile(frac, src, bkg, ratio=ratio, $
                 range=range, err_Ex=err_Ex, $
                 Qx=Qx, err_Qx=err_Qx, $
                 nofixerror=nofixerror, $
                 QDx=QDx, QDy=QDy, $
                 err_QDx=err_QDx, err_QDy=err_QDy)

   print,'----------------------'
   print,'src   ',n_elements(src)
   print,'bkg   ',n_elements(bkg)
   print,'net   ',n_elements(src)-ratio*n_elements(bkg)
   print,'ratio ',ratio
   print,'frac  ',frac
   print,'Ex    ',Ex
   print,'err_Ex',err_Ex
   print,'Qx    ',Qx
   print,'err_Qx',err_Qx
   print,'QDx   ',QDx, err_QDx
   print,'QDy   ',QDy, err_QDy

   QDx_array_gt_25(j)=QDx
   QDy_array_gt_25(j)=QDy
   QDxerr_array_gt_25(j, *)=err_QDx
   QDyerr_array_gt_25(j, *)=err_QDy
   j = j+1
endfor
; plot the data point
oploterror, QDx_array_lt_25, QDy_array_lt_25, QDxerr_array_lt_25(*,0), QDyerr_array_lt_25(*,0), /nohat, ERRCOLOR='green', PSYM=3
oploterror, QDx_array_gt_25, QDy_array_gt_25, QDxerr_array_gt_25(*,0), QDyerr_array_gt_25(*,0), /nohat, ERRCOLOR='red', PSYM=3


device,/close


end
