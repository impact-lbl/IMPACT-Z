import sys
import re
from collections import defaultdict
from copy import deepcopy
from math import sqrt, pi
import math
import const
from lattice_parser import lattice_parser

nest_dict = lambda: defaultdict(nest_dict) # to define a['key1']['key2'] = value
                 

# class: impactz_parser
#======================
class impactz_parser(lattice_parser):
    '''
    transform lte.impz => ImpactZ.in.
    
    '''
    def __init__(self, fileName, lineName):
        lattice_parser.__init__(self,fileName, lineName)

        # get brief lte.impz lines
        self.lines = self.get_brieflines()
        
        self.control = {}
        self.beam = {}
        self.lattice = nest_dict()
        
        # initiate with default values
        self.__default_control()
        self.__default_beam()
        self.__default_lattice()
        
        # update with read in *.impz file
        self.control = self.update_control()
        self.beam    = self.update_beam()
        self.lattice = self.update_trackline()

        # some parameters of the beam
        self.Scxl = None
        self.gam0 = None
        self.gambet0 = None
        self.bet0 = None
        self.mass = None
        self.update_paras()
        
    def write_impactzin(self):
        '''
        generate ImpactZ.in file.
        '''
        control_file = self.impactzin_control()
        lattice_file = self.impactzin_lattice()
        
        with open('ImpactZ.in','w') as f:
            f.write(control_file)
            f.write(lattice_file)

        f.close()     

    # default control, beam, lattice dicts
    #==============================================================================
    def __default_control(self):    
        self.control['CORE_NUM_T'] = 1  #transverse direction computor core number      
        self.control['CORE_NUM_L'] = 1  #longitudinal direction 
        self.control['INTEGRATOR'] = 1  # 1 for linear map, 2 for lorentz integrator
        self.control['ERROR'] = 0
        self.control['OUTPUT_RATIO'] = 1
        self.control['MESHX'] = 64
        self.control['MESHY'] = 64
        self.control['MESHZ'] = 64
        self.control['FLAGBC']= 1
        self.control['X_PIPE_WIDTH']= 0.014
        self.control['Y_PIPE_WIDTH']= 0.014
        self.control['PERIOD_LEN']= 0.1
        self.control['RESTART']   = 0
        self.control['SUB_CYCLE'] = 0
        self.control['KINETIC_ENERGY'] = 0;  # kinetic energy W, E=W+E0
        self.control['FREQ_RF_SCALE'] = 2.856e9
        self.control['INI_PHASE'] = 0.0

        self.control['STEPS_PERMETER'] = 0
        self.control['MAPS_PERMETER'] = 0
        self.control['SAMPLE_OUT'] = 1e5
        self.control['SLICE_BIN'] = 128
        self.control['PIPE_RADIUS'] = 0.014

        # turn all para values to str data type
        for key in self.control:
            self.control[key] = str( self.control[key] )

    def __default_beam(self):
        self.beam['MASS'] = const.electron_mass
        self.beam['CHARGE'] = -1.0
        self.beam['DISTRIBUTION_TYPE'] = 2
        self.beam['NP']   = int(1e3)
        self.beam['TOTAL_CHARGE'] = 1.0e-9 #[C]
        
        self.beam['SIGX'] = 0.0     #[m]
        self.beam['SIGY'] = 0.0     #[m]      
        self.beam['SIGZ'] = 0.0     #[m]
        
        self.beam['SIGPX'] = 0.0    #sig_gambetx/gambet0 [rad]
        self.beam['SIGPY'] = 0.0    #sig_gambety/gambet0 [rad]
        self.beam['SIGE']  = 0.0    # (E-E0)  [eV] 
        
        # for twiss parameters settings
        self.beam['EMIT_X'] = 0.0
        self.beam['EMIT_NX'] = 0.0   
        self.beam['BETA_X'] = 1.0
        self.beam['ALPHA_X'] = 0.0

        self.beam['EMIT_Y'] = 0.0
        self.beam['EMIT_NY'] = 0.0   
        self.beam['BETA_Y'] = 1.0
        self.beam['ALPHA_Y'] = 0.0

        self.beam['EMIT_Z'] = 0.0   # deg MeV
        self.beam['BETA_Z'] = 1.0   # deg/MeV
        self.beam['ALPHA_Z'] = 0.0  # 1

        # turn all para values to str data type
        for key in self.beam:
            self.beam[key] = str(self.beam[key])

    def __default_lattice(self):
        # drift, 0 element
        self.lattice['DRIFT']['L'] = 0.0
        self.lattice['DRIFT']['STEPS'] = 0
        self.lattice['DRIFT']['MAPS'] = 0
        self.lattice['DRIFT']['PIPE_RADIUS'] = 0.0

        # quad, 1 element
        self.lattice['QUAD']['L'] = 0.0
        self.lattice['QUAD']['STEPS'] = 0
        self.lattice['QUAD']['MAPS'] = 0
        self.lattice['QUAD']['K1'] = 0.0 
        self.lattice['QUAD']['GRAD'] = 0.0 
        self.lattice['QUAD']['PIPE_RADIUS'] = 0.0
        self.lattice['QUAD']['DX'] = 0.0
        self.lattice['QUAD']['DY'] = 0.0
        self.lattice['QUAD']['ROTATE_X'] = 0.0 #[rad]
        self.lattice['QUAD']['ROTATE_Y'] = 0.0
        self.lattice['QUAD']['ROTATE_Z'] = 0.0

        # BEND, 4 element
        self.lattice['BEND']['L'] = 0.0
        self.lattice['BEND']['STEPS'] = 0
        self.lattice['BEND']['MAPS'] = 0
        self.lattice['BEND']['ANGLE'] = 0.0
        self.lattice['BEND']['E1'] = 0.0
        self.lattice['BEND']['E2'] = 0.0
        self.lattice['BEND']['K1'] = 0.0
        self.lattice['BEND']['PIPE_RADIUS'] = 0.0 #half gap
        self.lattice['BEND']['H1'] = 0.0
        self.lattice['BEND']['H2'] = 0.0
        self.lattice['BEND']['FINT'] = 0.5
        self.lattice['BEND']['DX'] = 0.0
        self.lattice['BEND']['DY'] = 0.0
        self.lattice['BEND']['ROTATE_X'] = 0.0
        self.lattice['BEND']['ROTATE_Y'] = 0.0
        self.lattice['BEND']['ROTATE_Z'] = 0.0
        self.lattice['BEND']['CSR'] = 0

        # RFCW 
        self.lattice['RFCW']['L'] = 0.0
        self.lattice['RFCW']['STEPS'] = 0
        self.lattice['RFCW']['MAPS'] = 0
        self.lattice['RFCW']['VOLT'] = 0.0
        self.lattice['RFCW']['GRADIENT'] = 0.0
        self.lattice['RFCW']['PHASE'] = 0.0
        self.lattice['RFCW']['FREQ'] = 2.856e9
        self.lattice['RFCW']['PIPE_RADIUS'] = 0.0
        self.lattice['RFCW']['DX'] = 0.0
        self.lattice['RFCW']['DY'] = 0.0
        self.lattice['RFCW']['ROTATE_X'] = 0.0
        self.lattice['RFCW']['ROTATE_Y'] = 0.0
        self.lattice['RFCW']['ROTATE_Z'] = 0.0
        self.lattice['RFCW']['WAKEFLAG'] = -1
        self.lattice['RFCW']['WAKEFILE_ID'] = None

        # 101 DTL element 
        self.lattice['DTL']['L'] = 0.0
        self.lattice['DTL']['STEPS'] = 0
        self.lattice['DTL']['MAPS'] = 0
        self.lattice['DTL']['SCALE'] = 1.0
        self.lattice['DTL']['FREQ']  = 324e6
        self.lattice['DTL']['PHASE'] = 0.0     #driven phase, deg in sin  
        self.lattice['DTL']['ID']  = 100
        self.lattice['DTL']['PIPE_RADIUS'] = 0.0        
        self.lattice['DTL']['LQ1'] = 0.0    
        self.lattice['DTL']['GRAD1'] = 0.0        
        self.lattice['DTL']['LQ2'] = 0.0        
        self.lattice['DTL']['GRAD2'] = 0.0        
        self.lattice['DTL']['DX_Q'] = 0.0
        self.lattice['DTL']['DY_Q'] = 0.0
        self.lattice['DTL']['ROTATE_X_Q'] = 0.0
        self.lattice['DTL']['ROTATE_Y_Q'] = 0.0
        self.lattice['DTL']['ROTATE_Z_Q'] = 0.0
        self.lattice['DTL']['DX_RF'] = 0.0
        self.lattice['DTL']['DY_RF'] = 0.0
        self.lattice['DTL']['ROTATE_X_RF'] = 0.0
        self.lattice['DTL']['ROTATE_Y_RF'] = 0.0
        self.lattice['DTL']['ROTATE_Z_RF'] = 0.0

        # 104 SC element 
        self.lattice['SC']['L'] = 0.0
        self.lattice['SC']['STEPS'] = 0
        self.lattice['SC']['MAPS'] = 0
        self.lattice['SC']['SCALE'] = 1.0
        self.lattice['SC']['FREQ']  = 324e6
        self.lattice['SC']['PHASE'] = 0.0     #driven phase, deg in sin  
        self.lattice['SC']['ID']  = 100
        self.lattice['SC']['PIPE_RADIUS'] = 0.0      
        self.lattice['SC']['DX'] = 0.0
        self.lattice['SC']['DY'] = 0.0
        self.lattice['SC']['ROTATE_X'] = 0.0
        self.lattice['SC']['ROTATE_Y'] = 0.0
        self.lattice['SC']['ROTATE_Z'] = 0.0

        # 110 for user defined element: FIELDMAP
        self.lattice['FIELDMAP']['L'] = 0.0
        self.lattice['FIELDMAP']['STEPS'] = 0
        self.lattice['FIELDMAP']['MAPS'] = 0
        self.lattice['FIELDMAP']['SCALE'] = 1.0
        self.lattice['FIELDMAP']['FREQ']  = 324e6
        self.lattice['FIELDMAP']['PHASE'] = 0.0     #driven phase, deg in sin  
        self.lattice['FIELDMAP']['ID']  = 100
        self.lattice['FIELDMAP']['XRADIUS'] = 0.0        
        self.lattice['FIELDMAP']['YRADIUS'] = 0.0        
        self.lattice['FIELDMAP']['DX'] = 0.0
        self.lattice['FIELDMAP']['DY'] = 0.0
        self.lattice['FIELDMAP']['ROTATE_X'] = 0.0
        self.lattice['FIELDMAP']['ROTATE_Y'] = 0.0
        self.lattice['FIELDMAP']['ROTATE_Z'] = 0.0
        self.lattice['FIELDMAP']['DATATYPE'] = 1
        self.lattice['FIELDMAP']['COORDINATE'] = 2
        
        # -1 or -19, shift the centroid of beam to the axis origin point
        self.lattice['SHIFTCENTER']['OPTION'] = "ZDE"

        # watch, -2 and -8 element
        self.lattice['WATCH']['FILENAME_ID'] = 9999
        self.lattice['WATCH']['SAMPLE_FREQ'] = 0
        self.lattice['WATCH']['SLICE_BIN']   = 0
        self.lattice['WATCH']['SLICE_INFO'] = 1  
        self.lattice['WATCH']['COORD_INFO'] = 1  

        # -13, rectangular aperture 
        self.lattice['RCOL']['X_MAX'] =  0.04
        self.lattice['RCOL']['Y_MAX'] =  0.04
        self.lattice['RCOL']['X_MIN'] =  'None'
        self.lattice['RCOL']['Y_MIN'] =  'None'

        # -14, elliptical aperture 
        self.lattice['ECOL']['X_MAX'] = 0.04
        self.lattice['ECOL']['Y_MAX'] = 0.04

        # -18, ROTATE element
        self.lattice['ROTATE']['ANGLE']=0.0

        # scatter element, -20 element
        self.lattice['SCATTER']['DE']=0.0

        #turn all lattice elem values to string data type
        for elem in self.lattice.keys():
            for key in self.lattice[elem].keys():
                self.lattice[elem][key] = str(self.lattice[elem][key])

    # update the default control, beam, lattice with read-in lte.impz
    #==============================================================================    
    def update_control(self):
        '''
        update control:dict with read in para values.
        '''
        control_sec   = self.get_control_section()
        
        control = deepcopy( self.control )
        for key in control_sec.keys():
            if key in self.control.keys():               
                # update with read in values
                control[key] = control_sec[key]
            else:
                print('Unknown control item:',key,'=',control_sec[key])    
                sys.exit()
        return control
               
    def update_beam(self):
        '''
        update beam:dict with read in para values.
        '''
        beam_sec   = self.get_beam_section()
        
        beam = deepcopy( self.beam )
        for key in beam_sec.keys():
            if key in self.beam.keys():               
                # update with read in values
                beam[key] = beam_sec[key]
            else:
                print('Unknow beam item:',key,'=',beam_sec[key])    
                sys.exit()
        return beam

    def update_trackline(self):
        '''
        update the default track_line lattice para values with lte.impz
        '''
        trackline = self.get_lattice_section() # get the tracking line
        
        # MAP ELEGANT LATTICE TO lte.impz
        trackline = self.elegant2impz_lattice(trackline)
        
        j = 0
        for elem in trackline:
            # check if the element type is in self.lattice.keys, i.e. whether in
            # dict_keys(['DRIFT', 'QUAD', 'BEND', 'RFCW', 'WATCH'])          
                         
            # update the not-yet-setting lattice element parameters with the default 
            # values.
            if elem['TYPE'] in self.lattice.keys():
                tmp = deepcopy(self.lattice[elem['TYPE']])
                
                #NAME and TYPR for element in lte.impz are not in the self.lattice[elem['TYPE']].keys()
                table =  list(self.lattice[elem['TYPE']].keys())
                table.append('NAME')
                table.append('TYPE')
                for elem_para in elem.keys():
                    if elem_para not in table:
                        print("PAY ATTENTION: unknown element parameter",elem_para,"for",elem['NAME'],":",elem['TYPE'])
                        print("PROGRAM CONTINUE!")

                    tmp[elem_para] = elem[elem_para] 
                # replace trackline
                trackline[j] = tmp
                j += 1        
            else:
                print("Unknown element type in lattice section:",elem['TYPE'])
                sys.exit()       
        # turn all elem value to string data type    
        return trackline

    def update_paras(self):  
        Scxl = const.c_light/(2*math.pi*float(self.control['FREQ_RF_SCALE'])) 
        gam0 = (float(self.control['KINETIC_ENERGY'])+float(self.beam['MASS']))/float(self.beam['MASS'])
        gambet0 = sqrt(gam0**2-1.0)
        bet0 = gambet0/gam0

        #update self
        self.Scxl = Scxl
        self.gam0 = gam0
        self.gambet0 = gambet0
        self.bet0 = bet0
        self.mass = float(self.beam['MASS'])
     
    #write in ImpactZ.in    
    #==================================================================================================== 
    def impactzin_control(self):
        
        '''
        ImpactZ.in's control section
        '''
        # control section
        #----------------
        control_lines = [' !=================== \n']
        control_lines.append('! control section \n')
        control_lines.append('!=================== \n')
        
        # line-1 
        control_lines.append(self.control['CORE_NUM_T'])
        control_lines.append(self.control['CORE_NUM_L'])
        control_lines.append( '\n' )
        
        # line-2
        control_lines.append( '6' )
        Np = int(float(self.beam['NP']))
        control_lines.append( str(Np) )
        control_lines.append( self.control['INTEGRATOR'] )
        control_lines.append( self.control['ERROR'] )
        control_lines.append( self.control['OUTPUT_RATIO'] )
        control_lines.append( '\n' )
        
        # line-3
        control_lines.append( self.control['MESHX'] )
        control_lines.append( self.control['MESHY'] )
        control_lines.append( self.control['MESHZ'] )
        control_lines.append( self.control['FLAGBC'])
        control_lines.append(self.control['X_PIPE_WIDTH'])
        control_lines.append(self.control['Y_PIPE_WIDTH'])
        control_lines.append(self.control['PERIOD_LEN'])
        control_lines.append('\n')
        
        # line-4
        control_lines.append( self.beam['DISTRIBUTION_TYPE'] )
        control_lines.append( self.control['RESTART'] )
        control_lines.append( self.control['SUB_CYCLE'] )
        control_lines.append( '1 \n' )
        
        # line-5
        control_lines.append( str(Np) )
        control_lines.append( '\n' )
        
        # line-6
        # current averaged over FREQ_RF_SCALE
        current = float(self.beam['TOTAL_CHARGE'])*float(self.control['FREQ_RF_SCALE'])
        current = abs(current)  # current should be positive value even for electron beams,
                                # otherwise, negative current will turn off space charge.
        control_lines.append( str(current) )
        control_lines.append('\n')
        
        # line-7
        charge = float(self.beam['CHARGE'])
        q_m = 1/float(self.beam['MASS'])*charge  # q/m 
        control_lines.append(str(q_m))
        control_lines.append('\n')
        
        # line-8 to line-10
        # beam distribution  
        Scxl    = self.Scxl
        gam0    = self.gam0
        gambet0 = self.gambet0
        bet0    = self.bet0
        
        # (sigi,sigj,sigij) is different from what we know about RMS values,
        # only when sigij=0, sigi refers to RMS size. 
        # set sigij by default 0, for non-zero cases, use twiss-para.
        sigxxp=0.0
        sigyyp=0.0
        sigzdE=0.0

        # If twiss para is given, map twiss to sig
        # ----------------------------------------
        # X-PX
        emitx = float(self.beam['EMIT_X'])
        emit_nx = float(self.beam['EMIT_NX']) 
        if emit_nx != 0.0:
            emitx = emit_nx/gambet0
            
        if emitx != 0.0:
            betax = float(self.beam['BETA_X'])
            alphax = float(self.beam['ALPHA_X'])
            
            # attention: this is not RMS "sigx"
            sigx  = sqrt(emitx*betax/(1+alphax**2))
            sigxp = sqrt(emitx/betax)
            sigxxp= alphax/sqrt(1+alphax**2)
            
            # replace the default values
            self.beam['SIGX'] = str(sigx)
            self.beam['SIGPX'] = str(sigxp)
        
        # Y-PY
        emity = float(self.beam['EMIT_Y'])
        emit_ny = float(self.beam['EMIT_NY']) 
        if emit_ny != 0.0:
            emity = emit_ny/gambet0
            
        if emity != 0.0:           
            betay = float(self.beam['BETA_Y'])
            alphay = float(self.beam['ALPHA_Y'])
            
            sigy   = sqrt(emity*betay/(1+alphay**2))
            sigyp  = sqrt(emity/betay)
            sigyyp = alphay/sqrt(1+alphay**2)
             
            self.beam['SIGY']  = str(sigy)
            self.beam['SIGPY'] = str(sigyp)

        # Z-dE
        # deg-MeV
        emitz = float(self.beam['EMIT_Z'])

        if emitz != 0.0:           
            betaz = float(self.beam['BETA_Z'])
            alphaz = float(self.beam['ALPHA_Z'])
            
            sigz   = sqrt(emitz*betaz/(1+alphaz**2))
            sigdE  = sqrt(emitz/betaz)
            sigzdE = alphaz/sqrt(1+alphaz**2)
            
            # From (deg, MeV) change to (m, eV)
            # phi=wt,dW=Ei-E0 [MeV] =>
            # z=-bet0 ct, dE=Ei-E0 [eV]
            sigz   = sigz/180*pi *Scxl*bet0
            sigdE  = sigdE*1e6
            sigzdE = -sigzdE  #PAY ATTENTION to the minus sign

            self.beam['SIGZ'] = str(sigz)
            self.beam['SIGE'] = str(sigdE)
 
        # change to IMPACT-Z coordinates definition
        sigX    = float(self.beam['SIGX'])/Scxl  
        sigPx   = float(self.beam['SIGPX'])*gambet0
        sigY    = float(self.beam['SIGY'])/Scxl  
        sigPy   = float(self.beam['SIGPY'])*gambet0        
        sigT    = float(self.beam['SIGZ'])/Scxl/bet0   # m2rad
        sigPt = float(self.beam['SIGE'])/float(self.beam['MASS']) #eV 2 1 
        # sigxxp, sigyyp, sigzdE, no unit para 
        sigTPt = sigzdE
        
        control_lines.append( str(sigX) )
        control_lines.append( str(sigPx) )
        control_lines.append( str(sigxxp) )
        control_lines.append( '1.0 1.0 0.0 0.0 \n' )

        control_lines.append( str(sigY) )
        control_lines.append( str(sigPy) )
        control_lines.append( str(sigyyp) )
        control_lines.append( '1.0 1.0 0.0 0.0 \n' )
       
        control_lines.append( str(sigT)    )
        control_lines.append( str(sigPt) ) 
        control_lines.append( str(sigTPt)  ) 
        control_lines.append( '1.0 1.0 0.0 0.0 \n' )
        
        # line-11
        control_lines.append(str(current))
        control_lines.append( self.control['KINETIC_ENERGY'] )
        control_lines.append( self.beam['MASS'] )
        control_lines.append( self.beam['CHARGE'] )
        control_lines.append( self.control['FREQ_RF_SCALE'] )
        control_lines.append( self.control['INI_PHASE'] )
        control_lines.append( '\n' )

        control_lines = ' '.join(control_lines)
        return control_lines

    def impactzin_lattice(self):
        '''
        lattice section in ImpactZ.in file.        
        '''
        lte_lines = [' !=================== \n']
        lte_lines.append('! lattice lines \n')
        lte_lines.append('!=================== \n')        
               
        for elem in self.lattice:
            if elem['TYPE'] == 'DRIFT':
                if float(elem['L']) == 0:
                    # ignore 0.0 length drift
                    pass
                else:
                    self._set_steps_maps_radius(elem)

                    lte_lines.append( elem['L'] )
                    lte_lines.append( elem['STEPS'] )
                    lte_lines.append( elem['MAPS'] )
                    lte_lines.append('0')
                    lte_lines.append(elem['PIPE_RADIUS'])
                    lte_lines.append('/ \n')

            elif elem['TYPE'] == 'QUAD':
                self._set_steps_maps_radius(elem)

                lte_lines.append(elem['L'])
                lte_lines.append(elem['STEPS'])
                lte_lines.append(elem['MAPS'])
                lte_lines.append('1')
                lte_lines.append( elem['K1'] )
                lte_lines.append( '-1' )
                lte_lines.append(elem['PIPE_RADIUS'])
                lte_lines.append(elem['DX'])
                lte_lines.append(elem['DY'])
                lte_lines.append(elem['ROTATE_X'])
                lte_lines.append(elem['ROTATE_Y'])
                lte_lines.append(elem['ROTATE_Z'])
                lte_lines.append('/ \n')

            elif elem['TYPE'] == 'BEND':
                map_flag = self._get_bendmap_flag(elem)
                self._set_steps_maps_radius(elem)

                lte_lines.append(elem['L'])
                lte_lines.append(elem['STEPS'])
                lte_lines.append(elem['MAPS'])
                lte_lines.append('4')
                lte_lines.append(elem['ANGLE'])
                lte_lines.append(elem['K1']) 
                lte_lines.append( map_flag )
                lte_lines.append(elem['PIPE_RADIUS'])
                lte_lines.append(elem['E1']) 
                lte_lines.append(elem['E2']) 
                lte_lines.append(elem['H1'])
                lte_lines.append(elem['H2'])
                lte_lines.append(elem['FINT'])
                lte_lines.append(elem['DX'])
                lte_lines.append(elem['DY'])
                lte_lines.append(elem['ROTATE_X'])
                lte_lines.append(elem['ROTATE_Y'])
                lte_lines.append(elem['ROTATE_Z'])
                lte_lines.append('/ \n')
                
            elif elem['TYPE'] == 'RFCW':
                if float(elem['VOLT'])==0:
                    gradient = 0.0
                else:
                    gradient = float(elem['VOLT'])/float(elem['L']) #V/m
                phase = float(elem['PHASE']) - 90 #cos for impactz, elegant sin function
                
                # add -41 wakefield element
                # ------------------------
                if elem['WAKEFLAG']=='-1':
                    pass
                else:
                    # read-in wakefile
                    lte_lines.append('0 0 1 -41 1.0')
                    lte_lines.append(elem['WAKEFILE_ID'])
                    lte_lines.append(elem['WAKEFLAG'])
                    lte_lines.append('/ \n')

                self._set_steps_maps_radius(elem)
                # cavity lines
                # ------------------------
                lte_lines.append(elem['L'])
                lte_lines.append(elem['STEPS'])
                lte_lines.append(elem['MAPS'])
                lte_lines.append('103')
                lte_lines.append(str(gradient))
                lte_lines.append(elem['FREQ'])
                lte_lines.append( str(phase) )
                lte_lines.append( '-1.0' )
                lte_lines.append(elem['PIPE_RADIUS'] )
                lte_lines.append('/ \n')
                       
                # wakefield stops at the exit of cavity
                # ------------------------
                if elem['WAKEFLAG']=='-1':
                    pass
                else:
                    # OFF wake after RFCW element
                    lte_lines.append('0 0 1 -41 1.0 41 -1 / \n')           

            elif elem['TYPE'] == 'DTL':
                self._set_steps_maps_radius(elem)
                phase = float(elem['PHASE']) - 90

                lte_lines.append(elem['L'])
                lte_lines.append(elem['STEPS'])
                lte_lines.append(elem['MAPS'])
                lte_lines.append('101')
                lte_lines.append(elem['SCALE'])
                lte_lines.append(elem['FREQ'])
                lte_lines.append( str(phase) )
                lte_lines.append( elem['ID'] )
                lte_lines.append(elem['PIPE_RADIUS'] )
              
                lte_lines.append( elem['LQ1'] )
                lte_lines.append( elem['GRAD1'] )
                lte_lines.append( elem['LQ2'] )
                lte_lines.append( elem['GRAD2'] )

                lte_lines.append(elem['DX_Q'])
                lte_lines.append(elem['DY_Q'])
                lte_lines.append(elem['ROTATE_X_Q'])
                lte_lines.append(elem['ROTATE_Y_Q'])
                lte_lines.append(elem['ROTATE_Z_Q'])
                lte_lines.append(elem['DX_RF'])
                lte_lines.append(elem['DY_RF'])
                lte_lines.append(elem['ROTATE_X_RF'])
                lte_lines.append(elem['ROTATE_Y_RF'])
                lte_lines.append(elem['ROTATE_Z_RF'])
                lte_lines.append('/ \n')

            elif elem['TYPE'] == 'SC':
                self._set_steps_maps_radius(elem)
                phase = float(elem['PHASE']) - 90

                lte_lines.append(elem['L'])
                lte_lines.append(elem['STEPS'])
                lte_lines.append(elem['MAPS'])
                lte_lines.append('104')
                lte_lines.append(elem['SCALE'])
                lte_lines.append(elem['FREQ'])
                lte_lines.append( str(phase) )
                lte_lines.append( elem['ID'] )
                lte_lines.append(elem['PIPE_RADIUS'] )
                lte_lines.append(elem['DX'])
                lte_lines.append(elem['DY'])
                lte_lines.append(elem['ROTATE_X'])
                lte_lines.append(elem['ROTATE_Y'])
                lte_lines.append(elem['ROTATE_Z'])
                lte_lines.append('/ \n')

            elif elem['TYPE'] == 'FIELDMAP':
                self._set_steps_maps_radius(elem)
                phase = float(elem['PHASE']) - 90
                 
                if elem['XRADIUS']=='0.0':
                    elem['XRADIUS']=self.control['PIPE_RADIUS']
                if elem['YRADIUS']=='0.0':
                    elem['YRADIUS']=self.control['PIPE_RADIUS']

                lte_lines.append(elem['L'])
                lte_lines.append(elem['STEPS'])
                lte_lines.append(elem['MAPS'])
                lte_lines.append('110')
                lte_lines.append(elem['SCALE'])
                lte_lines.append(elem['FREQ'])
                lte_lines.append( str(phase) )
                lte_lines.append( elem['ID'] )
                lte_lines.append(elem['XRADIUS'] )
                lte_lines.append(elem['YRADIUS'] )
                lte_lines.append(elem['DX'])
                lte_lines.append(elem['DY'])
                lte_lines.append(elem['ROTATE_X'])
                lte_lines.append(elem['ROTATE_Y'])
                lte_lines.append(elem['ROTATE_Z'])
                lte_lines.append(elem['DATATYPE'])
                lte_lines.append(elem['COORDINATE'])
                lte_lines.append('/ \n')

            elif elem['TYPE'] == 'SHIFTCENTER':
                if elem['OPTION']=='ZDE':
                    lte_lines.append('0 0 0 -19')
                    lte_lines.append('/ \n')
                elif elem['OPTION']=="XY":
                    lte_lines.append('0 0 0 -1')
                    lte_lines.append('/ \n')
                else:
                    print("ERROR: Not available option for SHIFTCENTER is given:",elem['OPTION'])
                    sys.exit()
                   
            elif elem['TYPE'] == 'WATCH':
                # if SAMPLE_FREQ not set in the lattice line, replace with the value in control section
                if elem['SAMPLE_FREQ']=='0':
                    Np = float(self.beam['NP'])
                    sample_out = float(self.control['SAMPLE_OUT'])
                    if Np <= sample_out:
                        elem['SAMPLE_FREQ'] = '1'
                    elif Np > sample_out:
                        freq = math.ceil(Np/sample_out)
                        elem['SAMPLE_FREQ'] = str(freq)
                
                if elem['COORD_INFO']=='1':
                    lte_lines.append('0 0')
                    lte_lines.append(elem['FILENAME_ID'])
                    lte_lines.append('-2')
                    lte_lines.append('0.0')
                    lte_lines.append(elem['SAMPLE_FREQ'])
                    lte_lines.append('/ \n')
                
                # whether to add -8 element
                if elem['SLICE_BIN']=='0':
                    elem['SLICE_BIN']=self.control['SLICE_BIN']

                if elem['SLICE_INFO'] == '0':
                    pass
                
                elif elem['SLICE_INFO'] == '1':
                    lte_lines.append('0 0')
                    lte_lines.append( str(int(elem['FILENAME_ID']) +10000) )
                    lte_lines.append('-8')
                    lte_lines.append(elem['SLICE_BIN'])
                    lte_lines.append('/ \n')    
                else:
                    print('Unknown flag for SLICE_INFO, it should be 0 or 1.')
                    sys.exit()

            elif elem['TYPE']=='RCOL':
                if elem['X_MIN']=='None':
                    elem['X_MIN']=''.join(['-',elem['X_MAX']])
                if elem['Y_MIN']=='None':
                    elem['Y_MIN']=''.join(['-',elem['Y_MAX']])
                lte_lines.append('0 0 0 -13 0.014')
                lte_lines.append(elem['X_MIN'])
                lte_lines.append(elem['X_MAX'])
                lte_lines.append(elem['Y_MIN'])
                lte_lines.append(elem['Y_MAX'])
                lte_lines.append('/ \n')

            elif elem['TYPE']=='ECOL':
                lte_lines.append('0 0 0 -14 0.014')
                lte_lines.append(elem['X_MAX'])
                lte_lines.append(elem['Y_MAX'])
                lte_lines.append('/ \n')

            elif elem['TYPE']=='ROTATE':
                lte_lines.append('0 0 0 -18 0.014')
                lte_lines.append(elem['ANGLE'])
                lte_lines.append('/ \n')

            elif elem['TYPE'] == 'SCATTER':
                lte_lines.append('0 0 0 -20 0.014')
                lte_lines.append(elem['DE'])
                lte_lines.append('/ \n')


            else:
                print("NOT AVAILABLE ELEMENT TYPE:",elem['TYPE'])
                sys.exit()
       
        lte_lines = ' '.join(lte_lines)
        return lte_lines    
   
    # sub-funcs 
    #===============================================================================
   
    def get_control_section(self):
        '''
        get control section.
        '''
        lines = deepcopy(self.lines)
        
        j=0
        for line in lines:
           lines[j]=line.replace(' ','') 
           j=j+1
           
        # get control section
        pattern1 = re.compile(r'^&control$',re.I)  #not case sensitive
        pattern2 = re.compile(r'^&end$',re.I)
        
        j1,j2 = self._get_index(pattern1, pattern2, lines)
        
        control_lines = lines[j1+1:j2]  # ignored 1st and last element, i.e. &control and &end
        
        control = {}
        for line in control_lines:
            tmp = re.split(';|,',line) # , and ; are both supported
        
            # remove white space
            while '' in tmp:
                tmp.remove('')
            
            for j in tmp:
                tmp2 = j.split('=')
                
                name = tmp2[0].upper()
                value = tmp2[1]
                
                # in case math expression, such as:
                # FREQ_RF_SCALE = 2.998e8/2/pi
                try:
                    eval(value)
                except:
                    pass
                else:
                    value = eval(value.lower())
                    value = str(value)  #back to string
                
                control[name] = value                
                
        return control         
        
    def get_beam_section(self):
        '''
        get beam section.
        '''
        #use deepcopy to keep self.lines unchanged, as white space should be kept in
        #rpn expressions for lattice section  
        lines = deepcopy(self.lines)
        j=0
        for line in lines:
            lines[j] = line.replace(' ','')
            j=j+1
            
        # get lattice section
        pattern1 = re.compile(r'^&beam$',re.I)  #not case sensitive
        pattern2 = re.compile(r'^&end$',re.I)
        
        j1,j2 = self._get_index(pattern1, pattern2, lines)
        
        beam_lines = lines[j1+1:j2]  # ignored 1st and last element, i.e. &beam and &end

        beam = {}
        for line in beam_lines:
            tmp = re.split(';|,',line) # , and ; are both supported
        
            # remove white space
            while '' in tmp:
                tmp.remove('')
            
            for j in tmp:
                tmp2 = j.split('=')
                
                name = tmp2[0].upper()
                value = tmp2[1]
                
                # in case math expression, such as:
                #    total_charge=20*5e-3/2.998e8
                try:
                    eval(value)
                except:
                    pass
                else:
                    value = eval(value.lower())
                    value = str(value)  #back to string
                
                beam[name] = value
                                         
        return beam  

    def get_lattice_section(self):
        '''
        get lattice section. 

        '''   
        lines = self.lines
        
        # get lattice section
        pattern1 = re.compile(r'^&lattice$',re.I)  #not case sensitive
        pattern2 = re.compile(r'^&end$',re.I)
        
        j1,j2 = self._get_index(pattern1, pattern2, lines)
        
        lattice = lines[j1+1:j2]  # ignored 1st and last element, i.e. &lattice and &end
        
        # get the tracked line
        trackline = self.get_trackline(lattice)
        return trackline
               
    def _get_index(self, pattern1, pattern2, lines):
        '''
        get start and finished index of (&control,...,&end)    
        '''       
        j1 = int(1e6)
        j2 = int(1e6)
        cnt = 0
        for line in lines:
            if re.match(pattern1, line):  
                j1 = cnt
                break
            cnt += 1
                    
        cnt = j1
        for j in range(j1, len(lines)):
            if re.match(pattern2,lines[j]):
                j2 = cnt  
                break                  
            cnt += 1    
            
        if j1 == int(1e6):
            print(pattern1,"not found. Input file is wrong.")
            sys.exit()
        elif j2==int(1e6):
            print(pattern2,"not found. Input file is wrong.")
            sys.exit()
                      
        return j1,j2


    def _set_steps_maps_radius(self, elem :dict):
        '''
        set maps and steps based on control settings
        '''
        steps = float(self.control['STEPS_PERMETER'])
        maps  = float(self.control['MAPS_PERMETER'])
        length = float(elem['L'])

        if elem['STEPS'] == '0':
            elem['STEPS'] = str(math.ceil(steps*length))
        else:
            steps = float(elem['STEPS'])

        if elem['MAPS']=='0':
            elem['MAPS'] = str(math.ceil(maps*length))
        else:
            maps = elem['MAPS']

        # in case element length is 0.0, steps cannot be 0
        # if steps=0, then elem['STEPS']='1'
        if length==0.0 or steps==0.0:
            elem['STEPS'] = '1'
        if length==0.0 or maps==0.0:
            elem['MAPS'] = '1'

        if elem['PIPE_RADIUS'] == '0.0' :
            elem['PIPE_RADIUS'] = self.control['PIPE_RADIUS']    

    def _get_bendmap_flag(self, elem :dict):

        if elem['CSR']=='0':
            flag = '150'
        else:
            flag = '250'
        return flag

    def elegant2impz_lattice(self,trackline):
        '''
        For elegant lattice, there are different definition or notations for
        accelerator elements, map these elements to lte.impz CONVENTION

        Returns
        -------
        trackline

        '''
        # map ELEGANT elements
        j=0  
        for elem in trackline:
            # map sbend to bend
            if elem['TYPE'] in ['CSRCSBEN', 'CSRCSBEND','SBEN','CSBEND']:
                elem['TYPE'] = 'BEND'

            # map quadrupole to quad 
            elif elem['TYPE'] in ['QUADRUPOLE']:
                elem['TYPE'] = 'QUAD'

            # map monitor, hkicker, vkicker, sextupole to drift, temporary 
            # length is kept
            elif elem['TYPE'] in ['DRIF', 'CSRDRIFT', \
                                'MONITOR','MONI','MARK', \
                                'KICK','HKICKER','VKICKER','HKICK','VKICK', \
                                'MALIGN','CLEAN', \
                                'SEXTUPOLE', \
                                'CHARGE']:
                elem['TYPE'] = 'DRIFT'
            
            trackline[j]['TYPE']=elem['TYPE']
            j=j+1
        return trackline
        
if __name__=='__main__':
        
    # debug
    # ======
    file_name = './debugExample/lte.impz'   
    line_name = 'line'
    
    lte = impactz_parser(file_name,line_name)
    lte.write_impactzin()    
    
        
    
