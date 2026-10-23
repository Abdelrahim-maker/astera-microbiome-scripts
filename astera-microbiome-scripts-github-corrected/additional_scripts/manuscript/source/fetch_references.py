import urllib.request,json,concurrent.futures,pathlib
refs={
'Fierer2017':'10.1038/nrmicro.2017.87','Six2000':'10.1016/S0038-0717(00)00179-6','Sun2018':'10.3389/fmicb.2018.00699','Schmidt2018':'10.1371/journal.pone.0192953','Sommermann2018':'10.1371/journal.pone.0195345','Hartman2018':'10.1186/s40168-017-0389-9','Bowles2017':'10.1111/1365-2664.12815','Jansa2002':'10.1007/s00572-002-0163-z','Kabir2005':'10.4141/P03-160','Cotrufo2013':'10.1111/gcb.12113','Lehmann2015':'10.1038/nature16069','Brewer2017':'10.1038/nmicrobiol.2016.198','Daims2015':'10.1038/nature16461','vanKessel2015':'10.1038/nature16459','Delgado2018':'10.1126/science.aap9516','McMurdie2013':'10.1371/journal.pone.0061217','Gloor2017':'10.3389/fmicb.2017.02224','Mallick2021':'10.1371/journal.pcbi.1009442','Anderson2001':'10.1111/j.1442-9993.2001.01070.pp.x','Anderson2006':'10.1111/j.1541-0420.2005.00440.x','AndersonWalsh2013':'10.1890/12-2010.1','BH1995':'10.1111/j.2517-6161.1995.tb02031.x','Rottjers2018':'10.1093/femsre/fuy030','Weiss2016':'10.1038/ismej.2015.235','Martino2019':'10.1128/mSystems.00016-19','Blondel2008':'10.1088/1742-5468/2008/10/P10008','Friedman2012':'10.1371/journal.pcbi.1002687','Kurtz2015':'10.1371/journal.pcbi.1004226','McMurdie2014':'10.1371/journal.pcbi.1003531','Warton2012':'10.1111/j.2041-210X.2011.00127.x',
}
out=pathlib.Path('manuscript/source');raw=pathlib.Path('manuscript/audit/reference_metadata');raw.mkdir(exist_ok=True)
def get(kv):
 k,doi=kv
 try:
  req=urllib.request.Request('https://api.crossref.org/works/'+urllib.parse.quote(doi,safe=''),headers={'User-Agent':'ManuscriptReferenceVerification/1.0'})
  data=json.load(urllib.request.urlopen(req,timeout=45))['message'];(raw/(k+'.json')).write_text(json.dumps(data,indent=2))
  year=data.get('published',data.get('issued'))['date-parts'][0][0]
  if k=='Brewer2017':year=2017
  c={'id':k,'type':'article-journal','title':data['title'][0],'author':[{'family':a.get('family',''),'given':a.get('given','')} for a in data.get('author',[])],'issued':{'date-parts':[[year]]},'container-title':data.get('container-title',[''])[0],'volume':data.get('volume',''),'issue':data.get('issue',''),'page':data.get('page',data.get('article-number','')),'DOI':doi,'URL':'https://doi.org/'+doi}
  print(k,c['title'],flush=True);return c
 except Exception as e:print('FAILED',k,str(e),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:cs=[c for c in ex.map(get,refs.items()) if c]
(out/'references.json').write_text(json.dumps(cs,indent=2,ensure_ascii=False))
