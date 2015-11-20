; IDL Version 5.4 (sunos sparc)
; From Hugh Hudson
; Date: Sun Jul 20 11:45:38 2003

function break_time, t_in, qdebug = qdebug
  if n_elements(t_in) eq 0 then t_in = !stime
  extime = anytim(t_in,/ext)
  year = strtrim(string(extime(6)),1)
  month = strtrim(string(extime(5)),1)
  if strlen(month) eq 1 then month = '0'+month
  day = strtrim(string(extime(4)),1)
  if strlen(day) eq 1 then day = '0'+day
  hour = strtrim(string(extime(0)),1)
  if strlen(hour) eq 1 then hour = '0'+hour
  minute = strtrim(string(extime(1)),1)
  if strlen(minute) eq 1 then minute = '0'+minute
  second = strtrim(string(extime(2)),1)
  if strlen(second) eq 1 then second = '0'+second
  filname = year+month+day+'_'+hour+minute+second
  if keyword_set(qdebug) then stop
  return, filname
end
