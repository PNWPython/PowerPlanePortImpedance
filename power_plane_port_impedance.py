
def xi(m,n):
    if m == 0 and n == 0:
        return 1
    elif m == 0 or n == 0:
        return 2**(1/2.0)
    else:
        return 2

def kay_lossless(relative_permittivity, freq_radian, speed_light = 299792458):
    '''
    SI Units
    '''
    kay = ( relative_permittivity ** ( 1 / 2.0 ) ) * freq_radian / speed_light
    
    return kay
   

def propagationConstant(relative_permittivity, freq_radian, speed_light = 299792458):
    '''
    SI Units
    '''
    kay = ( relative_permittivity ** ( 1 / 2.0 ) ) * freq_radian / speed_light
    
    return kay
    
def pi():
    import math
    return math.pi
    
def kaysubn(m,wx):
    
    k = ( m * pi() / wx )
    
    return k

def cosignFunc(port1, port2, width, ith, jth):
    import cmath
    
    output =  \
            cmath.cos( ith * pi() * port1[0] / width[0] ) *\
            cmath.cos( ith * pi() * port2[0] / width[0] ) *\
            cmath.cos( jth * pi() * port1[1] / width[1] ) *\
            cmath.cos( jth * pi() * port2[1] / width[1] )
            
    return output
    
def frequencySweep(start, stop, steps, type = 'dec'):
    import math
    arrayofpoints = []
    
    for point in range(0,steps + 1):
        if type == 'lin':
            arrayofpoints.append(float( float(start) + float(point) * ( float(stop) - float(start) ) / float(steps) ))
        elif type == 'dec':
            arrayofpoints.append(10.0 ** (math.log10(start) + float( float(point) * ( math.log10(stop) - math.log10(start) ) / float(steps) )))
    # print arrayofpoints
    return arrayofpoints

def losslessplanecalculation(ports, width, plane_height, frequencyBand, rel_permittivity):
    '''
    requires consistent units in input
    
    Assumes port dimensions are much smaller than smallest wavelength
    -->Requires checker, otherwise add sync terms to z(w)
    '''
    permeability = 4.0 * pi() * 1E-7 
    permittivity = 8.854187817E-12
    speed_light = 299792458
    
    # fmax = frequencyBand[-1]
    # print fmax
    # b = 1.0
    # a = 1.0
    # m_limit = int( fmax * 2.0 * b * rel_permittivity ** (1 / 2.0) / speed_light ) #Missing a factor
    # n_limit = int( fmax * 2.0 * a * rel_permittivity ** (1 / 2.0) / speed_light ) #Missing a factor
    m_limit = 20
    n_limit = 20
    # print m_limit, n_limit
    
    for frequency in frequencyBand:
        # print "\nFreq:{:2E}\t".format(frequency),
        print "\n{:2E}\t".format(frequency),
        rad_freq = frequency * 2 * pi()
        for port in range(0,len(ports)):
            for txport in range(port,len(ports)):
                # print "Port({},{})\t".format(port,txport),
                z_sum = 0.0 + 0.0j
                for m in range(0,m_limit + 1):
                    for n in range(0,n_limit + 1):
                        z_index = z_sum + cosignFunc(ports[port], ports[txport], width, m,n) * \
                        ( ( xi(m,n) ** 2.0 ) / ( width[0] * width[1] * ( kaysubn(m,width[0])**2.0 + kaysubn(n,width[1])**2.0  - kay_lossless( rel_permittivity, rad_freq )**2.0 ) ) )
                        z_sum = z_index
                
                # print "{:2E}\t".format(abs(z_sum)),
                print "{:2E}\t".format(abs(1j * rad_freq * permeability * plane_height * z_sum)),
    return

def complex_rel_permittivity(freq_hz, rel_permittivity, frequency_extraction, limits, loss_tangent):
    '''
    This needs to be plotted, and the limits determined 
    
    limits[1,0] determine the linear portion of the imaginary part of the constant. 
    '''
    
    import math, cmath
    
    # global limits
    
    def theta(freq):
        partial =  ( 10**(limits[1]) + 1.0j * freq ) / ( 10**(limits[0]) + 1.0j * freq )
        return complex(math.log(partial.real),math.copysign(math.log(abs(partial.imag)),partial.imag) )
        
    result = \
        ( 1 / ( ( limits[1] - limits[0] ) * math.log(10) ) ) * theta(freq_hz) *\
        -1 * loss_tangent * rel_permittivity * ( limits[1] - limits[0] ) * math.log(10) / ( theta(frequency_extraction).imag * 1.0j ) +\
        rel_permittivity * ( 1 + loss_tangent * ( theta(frequency_extraction).real / theta(frequency_extraction).imag * 1.0j ) )
    
    return result

def transmissionplanecalculation(ports, width, plane_height, frequencyBand, rel_permittivity, rel_permittivity_frequency_extraction, metal_conductivity, dielectric_conductivity, vdd_thickness, loss_tangent, rel_permeability = 1):
    '''
    requires consistent units in input
    
    Assumes port dimensions are much smaller than smallest wavelength
    -->Requires checker, otherwise add sync terms to z(w)
    
    Based on "Characterization of Power Distribution Networks", Novak, Miller
    '''
    import cmath, os
    
    permeability = 4.0 * pi() * 1E-7 
    permittivity = 8.854187817E-12
    speed_light = 299792458
    
    limits = (1,12)
    
    # fmax = frequencyBand[-1]
    # print fmax
    # b = 1.0
    # a = 1.0
    # m_limit = int( fmax * 2.0 * b * rel_permittivity ** (1 / 2.0) / speed_light ) #Missing a factor
    # n_limit = int( fmax * 2.0 * a * rel_permittivity ** (1 / 2.0) / speed_light ) #Missing a factor
    m_limit = 20
    n_limit = 20
    # print m_limit, n_limit
    prmtvyFile = open(os.path.join(os.path.dirname(os.path.realpath(__file__)),'permittivity_plot.csv'),'w')
    print prmtvyFile
    
    print "\n" 
    print "Freq\t\tPort00(R)\t\tPort00(I)\t\tPort10(R)\t\tPort10(I)\t\tPort11(R)\t\tPort11(I)"
    for frequency in frequencyBand:
        
        rad_freq = frequency * 2 * pi()

        skinDepth = ( 2 / ( rad_freq * metal_conductivity * rel_permeability * permeability ) )**(1 / 2.0)

        ksubm = (1.0 - 1.0j) / skinDepth
        
        zedsubs = ( 2 * ksubm / metal_conductivity ) / ( cmath.tan(ksubm * vdd_thickness) )

        zedOmega = zedsubs + 1.0j * rad_freq * rel_permeability * permeability * plane_height
        
        #This is not calculating the limits and loss_tangent needs to be adjusted as well
        complex_calculated_Er = complex_rel_permittivity(frequency, rel_permittivity, rel_permittivity_frequency_extraction, limits, loss_tangent)
   
        whyOmega = ( 1 / plane_height ) * ( dielectric_conductivity + 1.0j * rad_freq * permittivity * complex_calculated_Er )

        prmtvyFile.write("{:2E},{:2E},{:2E},\n".format(frequency,complex_calculated_Er.real,complex_calculated_Er.imag))

        
        print "\n{:2E}\t".format(frequency),
        for port in range(0,len(ports)):
            for txport in range(port,len(ports)):
                # print "Port({},{})\t".format(port,txport),
                z_sum = 0.0 + 0.0j
                for m in range(0,m_limit + 1):
                    for n in range(0,n_limit + 1):
                        z_index = z_sum + cosignFunc(ports[port], ports[txport], width, m,n) * \
                            ( (xi(m,n) ** 2.0 ) / ( width[0] * width[1] * ( kaysubn(m,width[0])**2.0 + kaysubn(n,width[1])**2.0 + zedOmega * whyOmega ) ) )
                            
                        z_sum = z_index
                        
                # print zedOmega * z_sum, "\t",
                print "{:2E}\t{:2E}\t".format((zedOmega * z_sum).real,(zedOmega * z_sum).imag),
                
    prmtvyFile.close()
    return
    
if __name__ == "__main__":

    #SI UNITS
    port_locations                         =  [(0,12.7e-2),(25.4e-2,12.7e-2)]     #meters
    planewidthXY                           =  (25.4e-2,25.4e-2)                   #meters
    plane_height                           =  0.79e-3                             #meters
    frequencyBand                          =  frequencySweep(1e2,1e11,300,'dec')  #Hertz
    rel_permittivity                       =  4.0
    rel_permittivity_frequency_extraction  =  1e9                                 #Hertz
    metal_conductivity                     =  5.8e7                               #Siemens?
    dielectric_conductivity                =  2e-24                               #Siemens?
    vdd_thickness                          =  35e-6                               #meters
    loss_tangent                           =  0.01
    rel_permeability                       =  1.0

    
    # Prints a table of Port impedances vs. Frequency
    # losslessplanecalculation([(0,0),],(25.4e-2,25.4e-2),0.79e-3,frequencySweep(1e2,1e11,10,'dec'),4.0) #Test Value Correct
    transmissionplanecalculation(port_locations, planewidthXY, plane_height, frequencyBand, rel_permittivity, rel_permittivity_frequency_extraction, metal_conductivity, dielectric_conductivity, vdd_thickness, loss_tangent, rel_permeability)



