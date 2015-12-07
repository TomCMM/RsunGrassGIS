
import grass.script as grass
import grass.script.array as garray
import pandas as pd
import datetime
import numpy as np

class IrradiationRsun():
    """
    DESCRIPTION
        Calculate and return the irradiation calculated with the model Rsun
    """
    def __init__(self, namerasterelev):
        self.env = grass.gisenv()
        self.namerasterelev = namerasterelev
        self.rasterelev = namerasterelev+ "@" + self.env.get('MAPSET')# complete name
 
    def display_region(self):
        r = grass.read_command("g.region", flags='p' )
        print r
 
    def display_env(self):
        for key, value in self.env.items(): 
            print key, value
 
    def slope_aspect(self):
        """
        DESCRIPTION
            Calculate the slope and the aspect
        INPUT
            ribeirao_elev: name of the elev file to use ex "ribeirao_elev"
            ribeirao_aspect: name of the aspect file to give ex "ribeirao_aspect"
            ribeirao_slope: name of the slope file to give ex "ribeirao_slope"
        RETURN
            raster file of slope and aspect
         
        """
        slope = self.namerasterelev +"_slope"
        aspect = self.namerasterelev +"_aspect"
         
        self.namerasterslope = slope
        self.namerasteraspect = aspect
         
        self.rasteraspect = aspect + "@" + self.env.get('MAPSET')# complete name
        self.rasterslope = slope + "@" + self.env.get('MAPSET')# complete name
 
        grass.run_command("r.slope.aspect", overwrite=True, elevation=self.rasterelev, aspect=aspect, slope=slope)
 
    def irradiance(self,time,day, alb=0.2, lin=3, flat = None, outpath=None):
        """
        DESCRIPTION
            Calculate the Irradiance
        INPUT
            Outpath: True, use the path to save the global irradiance array on a file
        """
 
        try:
            self.rasteraspect
            self.rasterslope
        except AttributeError:
            self.slope_aspect()
 
        if flat:
            rasteraspect = self.rasteraspectflat
            rasterslope = self.rasterslopeflat
        else:
            rasteraspect = self.rasteraspect
            rasterslope = self.rasterslope
             
 
 
        grass.run_command("r.sun", 
                            flags="s",
                            overwrite=True,
                            verbose=False,
                            alb=alb,
                            lin=lin,
                            elevin=self.rasterelev, 
                            aspin=rasteraspect, 
                            slopein=rasterslope,
                            time=time,
                            beam_rad="beam_rad",
                            diff_rad="diff_rad",
                            refl_rad="refl_rad",
                            glob_rad="glob_rad",
                             
                            day = day
                            )
         
        self.rasterbeam = "beam_rad"+ "@" + self.env.get('MAPSET')# complete name
         
         
        self.rasterdiff = "diff_rad"+ "@" + self.env.get('MAPSET')# complete name
        self.rasterrefl = "refl_rad"+ "@" + self.env.get('MAPSET')# complete name
        self.rasterglob = "glob_rad"+ "@" + self.env.get('MAPSET')# complete name
        
        if outpath:
            irrglob = garray.array()
            irrglob.read(self.rasterelev)
            filename = outpath + str(day) + "_" + str(time)
            np.savetxt(filename, irrglob, delimiter=',') 
            print "Saved at :" + outpath
            
 
    def get_values(self, mapname, east, north):
        """
        DESCRIPTION
            return value of a point of a given map
        EXAMPLE
            grass.run_command("r.what", flags="fn", input="ASTGTM2_S23W047_dem_crop@PERMANENT", east_north=[-46.256593,-22.878141])
        """
        output = grass.read_command("r.what", flags="fn", input=mapname, east_north=[east,north])
        output = output.split("|")
        value = output[-2]
        return value
 
    def _new_map_horizontal_slope(self, north=None, east=None):
        """
        DESCRIPTION
            Write a new map with an horizontal slope and aspect
        TODO
            compare the results with the non horizontal
        """
        # get the value of the position to be change
#         self.get_values()
 
        try:
            self.rasterslope
            self.rasteraspect
        except AttributeError:
            self.slope_aspect()
 
 
        roundvalue = 4
        value = np.float64(0.00000000000000000001) # can not set 0 ?????
         
        elev_point = np.array(float(self.get_values(self.rasterelev, north=north, east=east)))
        slope_point = np.array(float(self.get_values(self.rasterslope, north=north, east=east)))
        aspect_point = np.array(float(self.get_values(self.rasteraspect, north=north, east=east)))
 
        elev_point =np.round(elev_point,roundvalue)
        slope_point =np.round(slope_point,roundvalue)
        aspect_point =np.round(aspect_point,roundvalue)
 
        elev = garray.array()
        slope = garray.array()
        aspect = garray.array()
 
#         print self.namerasterelev
#         print self.namerasterslope
#         print self.namerasteraspect
#         
        elev.read(self.namerasterelev)
        slope.read(self.namerasterslope)
        aspect.read(self.namerasteraspect)
         
        arrayelev = np.round(elev,roundvalue)
        arrayslope = np.round(slope,roundvalue)
        arrayaspect = np.round(aspect,roundvalue)
 
#         print np.round(elev,roundvalue)
#         print elev_point
#          
#         print np.round(slope,roundvalue)
#         print slope_point
#          
#         print np.round(aspect,roundvalue)
#         print aspect_point
#  
#         print zip(*np.where(arrayelev == elev_point))
#         print zip(*np.where(arrayslope == slope_point))
#         print zip(*np.where(arrayaspect == aspect_point))
 
 
        i_elev = zip(*np.where(arrayelev == elev_point))
        i_slope = zip(*np.where(arrayslope == slope_point))
        i_aspect = zip(*np.where(arrayaspect == aspect_point))
 
        newset = set(i_elev).intersection(i_slope)
        newset = newset.intersection(i_aspect)
 
        index = list(newset)[0]
        index = [index[0], index[1]]
        print index
 
        slope[index[0], index[1]] = value
        aspect[index[0], index[1]] = value
 
#         
#         print "="*120
#         print "New valor set"
 
 
        outslope = self.namerasterslope+"_flat"
        outaspect = self.namerasteraspect+"_flat"
         
        self.namerasterslopeflat = outslope
        self.namerasteraspectflat = outaspect
         
        self.rasterslopeflat = self.namerasterslopeflat + "@" + self.env.get('MAPSET')# complete name
        self.rasteraspectflat = self.namerasteraspectflat + "@" + self.env.get('MAPSET')# complete name
 
        slope.write(outslope, overwrite = True)
        aspect.write(outaspect, overwrite = True)
#         
#         print self.get_values(self.rasterslopeflat, north=north, east=east)
#         
         
        print "WRITE THE NEW SLOPE AND ASPECT HORIZONTAL"
 
    def get_serie_irradiance(self, days, hours, lag=0,lin = 3, east=None, north=None, years = [2014,2015],outpath=None, flat=None):
        """
        DESCRIPTION
            return a pandas dataframe with the component of the radiation
            at a specified point
        INPUT
            years: just to put a year at the index datetime 
            then the data are easily comparable to the observations
            lag: apply a lag to the measure
            lin: Linketurbidity factor for the calcul of the divergence
            years: years to perform the simulation
        """

        if not isinstance(east,list):
            east = [east]
            north = [north]
        
        days_column = []
        hours_column = []
        df_names = ['beam_df', "diff_df", "refl_df", "glob_df", "slope_df", "aspect_df", "elev_df"]
        beam_df = pd.DataFrame()
        diff_df = pd.DataFrame()
        refl_df = pd.DataFrame()
        glob_df = pd.DataFrame()
        slope_df = pd.DataFrame()
        aspect_df = pd.DataFrame()
        elev_df = pd.DataFrame()

        rasterelev = self.rasterelev
        
        try:
            self.rasterslope
            self.rasteraspect
        except AttributeError:
            self.slope_aspect()

        if flat:
            print "="*120
            print "FLAT"
            print "="*120
            self._new_map_horizontal_slope(north=north,east=east)
            rasteraspect = self.rasteraspectflat
            rasterslope = self.rasterslopeflat
        else:
            rasteraspect = self.rasteraspect
            rasterslope = self.rasterslope
         
        for day in days:
            print day
            for hour in hours:
                time = np.array(hour+lag)
                self.irradiance(day=day, time=time, flat = flat, lin=lin) # Compute the irradiance
                days_column.append(day)
                hours_column.append(hour)

                beam = []
                diff = []
                refl = []
                glob = []
                slope = []
                aspect = []
                elev = []
                
                
                for e,n in zip(east, north):
                    beam.append(self.get_values(self.rasterbeam, east=e, north=n))
                    diff.append(self.get_values(self.rasterdiff, east=e, north=n))
                    refl.append(self.get_values(self.rasterrefl, east=e, north=n))
                    glob.append(self.get_values(self.rasterglob, east=e, north=n))
                    slope.append(self.get_values(rasterslope, east=e, north=n))
                    aspect.append(self.get_values(rasteraspect, east=e, north=n))
                    elev.append(self.get_values(rasterelev, east=e, north=n))


                beam_df = beam_df.append(pd.Series(beam), ignore_index=True)
                diff_df = diff_df.append(pd.Series(diff), ignore_index=True)
                refl_df = refl_df.append(pd.Series(refl), ignore_index=True)
                glob_df = glob_df.append(pd.Series(glob), ignore_index=True)
                slope_df = slope_df.append(pd.Series(slope), ignore_index=True)
                aspect_df = aspect_df.append(pd.Series(aspect), ignore_index=True)
                elev_df = elev_df.append(pd.Series(elev), ignore_index=True)

        dfs = [beam_df, diff_df, refl_df, glob_df, aspect_df, slope_df, elev_df]
        dfs_res = {}
        for df, name in zip(dfs, df_names):
            print df
            df['hours'] = hours_column
            df["days"] = days_column
            df['years'] = years[0]# THIS WILL NOT WORK FOR MORE THAN TWO YEARS BUT I AM TIRED
            for i in years[1:]:
                newdf = df.copy()
                newdf['years'] = i
                df = pd.concat([df, newdf])
              
            df['date'] = df.loc[:,['years', 'days','hours']].apply(self._dates, axis=1)
     
            df.index = df['date']
            del df['date']
            del df['years']
            del df['hours']
            del df['days']
            df.columns = [str(e)+str(n) for e,n in zip(east,north)]
            
            if outpath:
                df.to_csv("/home/thomas/Irradiance_rsun"+"_lin"+str(lin)+"_lag"+str(rsun_lag)+"_"+str(name)+".csv")
            else:
                dfs_res[name] = df

        if not outpath:
            return dfs_res
     
    def _dates(self,row):
        year = datetime.datetime(row['years'],1,1)
        timeday = datetime.timedelta(days=int(row['days']))
        hours = pd.DateOffset(hours=int(row['hours']))
 
        date = year + timeday + hours
        return date


if __name__ == '__main__':
    namerasterelev = "ASTGTM2_S23W047_dem_crop"

    irr = IrradiationRsun(namerasterelev)
    irr.display_region()
    irr.display_env()
    
    hours = np.arange(5,21)
    days = range(1,2)
    
    for day in days:
        for hour in hours:
            irr.irradiance(hour, day, 2, 2, flat=None, outpath="/home/thomas/")



#===============================================================================
# Get a serie of value
#===============================================================================
# # C4
# # C5
# # C6
# # C7
# # C8
# # C9
# # C10
# # C11
# # C12
# # C13
# # C14
# # C15
# # C16
# # C17
# # C18
# # C19
# 
#     north = [
#             -22.88097,
#             -22.88117,
#             -22.87792,
#             -22.87686,
#             -22.87411,
#             -22.87019,
#             -22.88331,
#             -22.88344,
#             -22.88628,
#             -22.88839,
#             -22.88914,
#             -22.88964,
#             -22.86303,
#             -22.86478,
#             -22.86414,
#             -22.86439
#             ]
#  
#     east = [
#             -46.249083,
#             -46.251667,
#             -46.252861,
#             -46.254528,
#             -46.256667,
#             -46.258833,
#             -46.246944,
#             -46.245861,
#             -46.243694,
#             -46.241278,
#             -46.238472,
#             -46.237139,
#             -46.24731,
#             -46.24394,
#             -46.23861,
#             -46.23614
#              ]
#     irr = IrradiationRsun(namerasterelev)
#     irr.display_region()
#     irr.display_env()
# 
#     r_sunlags = [-0.2]
#     lins = [2]
#     hours = np.arange(5,21)
#     days = range(1,355)
#     for rsun_lag in  r_sunlags:
#         for lin in lins:
#             df = irr.get_serie_irradiance(days, hours, east=east, north=north, lag = rsun_lag, lin=lin, outpath='/home/thomas')

