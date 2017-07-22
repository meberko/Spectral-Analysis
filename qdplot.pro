; JaeSub Hong, 2003-2005, version 1.7
; Please report any problem or suggestion at jaesub@head.cfa.harvard.edu
; 
; Calculate quantiles from a distribution
; Refer to 
;	J. Hong, E.M. Schlegel & J.E. Grindlay, 
;	2004 ApJ 614 p508 and references therein
; 
; QCCD diagram plotting routine
; require rdrdb.pro (included)

;========================================================================

pro qd_setwin, xr=xr, yr=yr, $
	xtit=xtit,ytit=ytit, $
	xty=xty,yty=yty, $
	xtick=xt, range=range, $
	stu=stu
; setting up window

if not keyword_set(xtit) then $
	xtit='log!d10!n(!s!a m !r!s!b1-m!r!s-!r !s-!r !s-!r -!3!n)' 
if not keyword_set(ytit) then $
	ytit='3 Q!d25!n/Q!d75!n'

plot,[0],[0],/nodata, $
 	xr=xr,xs=9,yr=yr,/ys,$
	xty=xty,yty=yty, $
	xtit=xtit,ytit=ytit

; labeled xtick
if not keyword_set(xt) then xt= ['0.6','1','2','3','4','5','6 keV','7'];

rl = range[1]-range[0]
m = (float(xt)-range[0])/rl
xv =alog10(m/(1.-m))
axis, xaxis=1,xticks=n_elements(xt)-1, xtickv=xv, xtickn=xt;,xticklen=0.01

; unlabeled small xtick

if keyword_set(stu) then begin
	min_m = (10.^xr[0]/(10.^xr[0]+1.))*rl+range[0]
	max_m = (10.^xr[1]/(10.^xr[1]+1.))*rl+range[0]

	min_m = fix(min_m/stu/10.-1)*stu*10
	max_m = fix(max_m/stu/10.+1)*stu*10
	n_sxt = (max_m-min_m)/stu
	if n_sxt gt 60. then begin
		print,'# of ticks should be less than 60'
		goto, done
	endif
	sxtv = findgen(n_sxt)*stu+min_m
	sxtv = (sxtv-range[0])/rl
	sxtv = alog10(sxtv/(1.-sxtv))
	sxt = strarr(n_sxt)+' '
	axis, xaxis=1,xticks=n_sxt-1, xtickv=sxtv, xtickn=sxt, xticklen=0.01
	done:
endif


end

function qd_readgrid, file_grid, comment=comment, n_data=n_data

; read grid
; if it's multiple, take the average
; assuming the identical format

n_data 	= rdrdb(file_grid[0],data=grid,comment=comment)
ntag 	= n_tags(grid)
n_file	= n_elements(file_grid)

for i=1, n_file-1 do begin
	n_data = rdrdb(file_grid[i],data=grid_,comment=comment)
	for j=3, ntag-1 do grid.(j) = grid.(j)+grid_.(j)
endfor

for j=3, ntag-1 do grid.(j) = grid.(j)/n_file

return, grid
end

pro qd_plotgrid, grid, comment=comment, $ 
	file_grid=file_grid, $ ; read the data
	range=range, $
	line=line, thick=thick, $
	col=colgrid, cid=cid

if not keyword_set(line) 	then line=0
if not keyword_set(colgrid) 	then colgrid=0
if not keyword_set(thick) 	then thick=!p.thick

if n_elements(thick) 	eq 1 	then thick=[thick,thick]
if n_elements(line) 	eq 1 	then line=[line,line]
if n_elements(colgrid) 	eq 1 	then colgrid=[colgrid,colgrid]

if keyword_set(file_grid) then $
	grid=qd_readgrid(file_grid, comment=comment)

if not keyword_set(range) then begin
	for i=0, n_elements(comment)-1 do begin
		pos=strpos(comment[i], '-range',0)
		if pos ge 0 then begin
			str_range = strtrim(strmid(comment[i],pos+7,20))
			goto,done
		endif
	endfor
	done:

	fields = strsplit(str_range,'[:,]',/extract,/reg)
	range = [float(fields[0]),float(fields[1])]
endif

rl = range[1]-range[0]

; this is safer, but keep compatiblity 
if keyword_set(cid) then begin
	Q25 = (grid.(cid[0])-range[0])/rl
	Q50 = (grid.(cid[1])-range[0])/rl
	Q75 = (grid.(cid[2])-range[0])/rl
endif else begin
	Q25 = (grid.E25-range[0])/rl
	Q50 = (grid.E50-range[0])/rl
	Q75 = (grid.E75-range[0])/rl
endelse

QDx = alog10(Q50/(1.-Q50))
QDy = 3.*Q25/Q75

gidt = strmid(grid.gridid,0,2)
gidn = fix(strmid(grid.gridid,2,3))

; now the first parameter
w=where(gidt eq 'gx',cw)
gidn1 = gidn[w]
QDx_ = QDx[w]
QDy_ = QDy[w]

for i=0, max(gidn1) do begin
	w=where(gidn1 eq i,cw)
	oplot, QDx_[w], QDy_[w], line=line[0], $
		col=colgrid[0],thick=thick[0];,psym=7
endfor

; next the second parameter
w=where(gidt eq 'gy',cw)
gidn2 = gidn[w]
QDx_ = QDx[w]
QDy_ = QDy[w]

for i=0, max(gidn2) do begin
	w=where(gidn2 eq i,cw)
	oplot, QDx_[w], QDy_[w], line=line[1], $
		col=colgrid[1],thick=thick[1];,psym=7
endfor

end

pro qd_plot_multigrid, file_grid, grid=grid, comment=comment, $
	type=type, col=col,$
	range=range, $
	line=line, thick=thick, cid=cid

if not keyword_set(type) then type='all'
if type eq 'all' then begin
	for i=0, n_elements(file_grid)-1 do begin
		grid=qd_readgrid(file_grid[i], comment=comment)
		qd_plotgrid, grid, comment=comment, col=col, $
			range=range, line=line, thick=thick, cid=cid
	endfor
endif else if type eq 'avg' then begin
	grid=qd_readgrid(file_grid, comment=comment)
	qd_plotgrid, grid, comment=comment, col=col, $
			range=range, line=line, thick=thick, cid=cid
endif

end

pro qd_labelgrid, par, grid, $
	label1=label1, label2=label2, $
	prefix1=prefix1, prefix2=prefix2, $
	suffix1=suffix1, suffix2=suffix2, $
	comment=comment, $ ; when grid data is given directly
	file_grid=file_grid, $ ; read the data
	range=range, $
	align=align, ori=ori, col=col, xo=xo, yo=yo, cid=cid

if keyword_set(file_grid) then $
	grid=qd_readgrid(file_grid, comment=comment)

if not keyword_set(align) then align=[0.0,0.0]
if not keyword_set(ori) then ori=[0.0,0.0]
if not keyword_set(xo) then xo=[0.0,0.0]
if not keyword_set(yo) then yo=[0.0,0.0]

if not keyword_set(range) then begin
	for i=0, n_elements(comment)-1 do begin
		pos=strpos(comment[i], '-range',0)
		if pos ge 0 then begin
			str_range = strtrim(strmid(comment[i],pos+7,20))
			goto,next
		endif
		pos=strpos(comment[i], '-par1',0)
		if pos ge 0 then begin
			str_par1 = strtrim(strmid(comment[i],pos+6,40))
			goto,next
		endif
		pos=strpos(comment[i], '-par2',0)
		if pos ge 0 then begin
			str_par2 = strtrim(strmid(comment[i],pos+6,40))
			goto,next
		endif
		next:

	endfor

	fields = strsplit(str_range,'[:,]',/extract,/reg)
	range = [float(fields[0]),float(fields[1])]

	if not keyword_set(prefix1) then prefix1=''
	if not keyword_set(prefix2) then prefix2=''
	if not keyword_set(suffix1) then suffix1=''
	if not keyword_set(suffix2) then suffix2=''

	fields = strsplit(str_par1,':',/extract)
	par1v = fields[1]
	par1v = prefix1+strsplit(par1v,',',/extract)+suffix1

	fields = strsplit(str_par2,':',/extract)
	par2v = fields[1]
	par2v = prefix2+strsplit(par2v,',',/extract)+suffix2

endif

rl = range[1]-range[0]

; this is safer, but keep compatiblity 
if keyword_set(cid) then begin
	Q25 = (grid.(cid[0])-range[0])/rl
	Q50 = (grid.(cid[1])-range[0])/rl
	Q75 = (grid.(cid[2])-range[0])/rl
endif else begin
	Q25 = (grid.E25-range[0])/rl
	Q50 = (grid.E50-range[0])/rl
	Q75 = (grid.E75-range[0])/rl
endelse

QDx = alog10(Q50/(1.-Q50))
QDy = 3.*Q25/Q75

gidt = strmid(grid.gridid,0,2)
gidn = fix(strmid(grid.gridid,2,3))

w=where(gidt eq 'gx',cw)
gidn1 = gidn[w]
grid1_ = grid[w].(1)
grid2_ = grid[w].(2)
QDx_ = QDx[w]
QDy_ = QDy[w]

for i=0, max(gidn1) do begin
	w=where(gidn1 eq i,cw)
	dum=min(abs(grid2_[w] - par[1]),ind)
	if n_elements(label1)-1 ge i then label_ = label1[i] $
	else label_ = par1v[i]
	xyouts, QDx_[w[ind]]+xo[0], QDy_[w[ind]]+yo[0], label_, $
		align=align[0], ori=ori[0],noclip=0,col=col[0]
endfor

w=where(gidt eq 'gy',cw)
gidn2 = gidn[w]
grid1_ = grid[w].(1)
grid2_ = grid[w].(2)
QDx_ = QDx[w]
QDy_ = QDy[w]

for i=0, max(gidn2) do begin
	w=where(gidn2 eq i,cw)
	dum=min(abs(grid1_[w] - par[0]),ind)
	if n_elements(label2)-1 ge i then label_ = label2[i] $
	else label_ = par2v[i]
	xyouts, QDx_[w[ind]]+xo[1], QDy_[w[ind]]+yo[1], label_, $
		align=align[1], ori=ori[1],noclip=0,col=col[1]
endfor


end

pro qdplot
	print, "qccd plot: version 1.7"
end

