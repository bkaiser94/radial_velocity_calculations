"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-09-10

Ideally this will contain the functions and methods that wave_cal.py needs to use for its extraction, but by placing 
the methods externally in here, I should guarantee uniform methods in my extractions and that I actually   
implement the things correctly...
This is probably not going to work.



"""
import numpy as np
import matplotlib.pyplot as plt



import spec_plot_tools as spt

def extract_spectrum(img_data, header, polynomials= [0,0], core_sides=2., bkg_core_sides= 10.,  bkg_shift=10., bkg_method= 'avg', bkg_poly_deg= 2., n_biases=0., trace_offset=0):
    band_inds= np.indices(img_data.shape)
    x_positions= band_inds[1,1]
    y_positions= band_inds[0]
    target_light= np.array([])
    target_noise2_list= np.array([])
    bkg_light= np.array([])
    bkg_noise2_list= np.array([])
    bkg_up_comb=np.array([])
    bkg_down_comb= np.array([])
    poly_curve_y = np.polyval(polynomials[0], x_positions)
    poly_curve_wavelength= np.polyval(polynomials[1], x_positions)
    upper_edges= np.polyval(polynomials[1], x_positions+0.5) #upper wavelength edges
    lower_edges= np.polyval(polynomials[1], x_positions-0.5) #lower_wavelength edges
    dlambda_vals= upper_edges-lower_edges
    
    def bkg_trace(x_positions, sign='minus'):
        #return  np.int_(poly_curve_y[x_positions]+bkg_shift)
        if sign=='minus':
            return  spt.discrete_int(poly_curve_y[x_positions])-bkg_shift
        elif sign=='plus':
            return spt.discrete_int(poly_curve_y[x_positions])+bkg_shift
    for x_pos in x_positions:
        #trace_vals=img_data[np.int_(poly_curve_y[x_pos]-core_sides):np.int_(poly_curve_y[x_pos]+core_sides+1),x_pos]
        trace_vals=img_data[spt.discrete_int(poly_curve_y[x_pos])-core_sides:spt.discrete_int(poly_curve_y[x_pos])+core_sides+1,x_pos]
        xsum= np.sum(trace_vals)
        #xsum= np.sum(img_data[np.int_(poly_curve_y[x_pos]-core_sides):np.int_(poly_curve_y[x_pos]+core_sides+1),x_pos]) #old way 2018-10-31
        target_light= np.append(target_light,[xsum])
       
        #up_bkg=img_data[np.int_(poly_curve_y[x_pos]+bkg_shift-bkg_core_sides):np.int_(poly_curve_y[x_pos]+bkg_shift+bkg_core_sides+1),x_pos]
        
        up_bkg=img_data[bkg_trace(x_pos, sign='plus')-bkg_core_sides:bkg_trace(x_pos, sign='plus')+bkg_core_sides+1,x_pos]
        #down_bkg= img_data[np.int_(poly_curve_y[x_pos]-bkg_shift-bkg_core_sides):np.int_(poly_curve_y[x_pos]-bkg_shift+bkg_core_sides+1),x_pos]
        down_bkg= img_data[bkg_trace(x_pos, sign='minus')-bkg_core_sides:bkg_trace(x_pos, sign='minus') +bkg_core_sides+1,x_pos]
        bkg_comb= np.append(up_bkg, down_bkg)
        if bkg_method=='poly':
            #up_bkg_coords= np.arange(bkg_trace(x_pos, sign='plus')-bkg_core_sides,bkg_trace(x_pos, sign='plus')+bkg_core_sides+1, 1)
            #down_bkg_coords= np.arange(bkg_trace(x_pos, sign='minus')-bkg_core_sides,bkg_trace(x_pos, sign='minus')+bkg_core_sides+1, 1)
            up_bkg_coords=y_positions[bkg_trace(x_pos, sign='plus')-bkg_core_sides:bkg_trace(x_pos, sign='plus')+bkg_core_sides+1,x_pos]
            down_bkg_coords= y_positions[bkg_trace(x_pos, sign='minus')-bkg_core_sides:bkg_trace(x_pos, sign='minus') +bkg_core_sides+1,x_pos]
            bkg_coords= np.append(up_bkg_coords, down_bkg_coords)
            
            trace_coords= np.arange(spt.discrete_int(poly_curve_y[x_pos])-core_sides, spt.discrete_int(poly_curve_y[x_pos])+core_sides+1, 1)
            #print('bkg_coords', bkg_coords)
            #print('bkg_comb', bkg_comb)
            bkg_polynomial= np.polyfit(bkg_coords, bkg_comb, bkg_poly_deg)
            if x_pos%1000 == 0:
                plt.plot(img_data[:,x_pos], label="Full cross section of flux", linestyle= 'none', marker='o')
                plt.plot(bkg_coords, bkg_comb, label='bkg points', linestyle='none', marker='o')
                rep_coords= np.linspace(np.nanmin(bkg_coords), np.nanmax(bkg_coords), 1000)
                plt.plot(rep_coords, np.polyval(bkg_polynomial, rep_coords), label='bkg fit')
                plt.plot(trace_coords, trace_vals, label='trace flux')
                plt.title("spext plot:"+  str(x_pos))
                plt.xlabel('y coord')
                plt.ylabel('counts')
                plt.legend()
                plt.show()
            bkg_calc_vals= np.polyval(bkg_polynomial, trace_coords)
            bkg_sum= np.sum(bkg_calc_vals)
        elif bkg_method== 'avg':
            bkg_sum= trace_vals.shape[0]*np.copy(np.mean(bkg_comb)) #take the mean of the bkg portion of the sky
            
            
        #plt.legend()
        #print "trace_vals.shape", trace_vals.shape
        #print "bkg_comb.shape", bkg_comb.shape
        target_noise2 = np.copy(xsum+trace_vals.shape[0]*header['RDNOISE']**2+trace_vals.shape[0]*header['RDNOISE']**2/n_biases)
        target_noise2_list= np.append(target_noise2_list, [target_noise2])
        bkg_noise2= trace_vals.shape[0]*np.copy(np.mean(bkg_comb)/bkg_comb.shape[0]+header['RDNOISE']**2/bkg_comb.shape[0]+header['RDNOISE']**2/(bkg_comb.shape[0]*n_biases))
        #bkg_sum= np.sum(img_data[np.int_(poly_curve_y[x_pos]+bkg_shift-core_sides):np.int_(poly_curve_y[x_pos]+bkg_shift+core_sides+1),x_pos])
        bkg_light= np.append(bkg_light,[bkg_sum])
        bkg_noise2_list= np.append(bkg_noise2_list, [bkg_noise2]) #list of noise values for a single pixel (resulting from the mean of the sky) for a given column
        bkg_up_comb= np.append(bkg_up_comb, [np.sum(up_bkg)])
        bkg_down_comb= np.append(bkg_down_comb, [np.sum(down_bkg)])
    #plt.plot(bkg_up_comb, label='bkg_up_comb')
    #plt.plot(bkg_down_comb, label='bkg_down_comb')
    #plt.plot(target_light, label='target_light')
    #plt.title('Before sky subtraction')
    #plt.legend()
    #plt.xlabel('pixel')
    #plt.ylabel('counts')
    #plt.show()
    #plt.plot(x_positions,target_light,'-')
    #plt.xlabel('x (pixel)')
    #plt.ylabel('Counts')
    #plt.title('Target Spectrum')
    #plt.show()
    #noise_spectrum = np.copy(np.sqrt(target_light + bkg_light + y_trace_width*header['RDNOISE'])) #old way 2018-10-31
    noise_spectrum= np.copy(np.sqrt(target_noise2_list+bkg_noise2_list)) #combination of noises of the background pixels and the target trace pixels.
    print "noise_spectrum.shape", noise_spectrum.shape
    target_light= target_light-bkg_light
    noise_spectrum = noise_spectrum/target_light #normalized noise values by the target spectrum, so now unitless.
    target_light= target_light/header['EXPTIME'] #converting to counts/s
    bkg_light= bkg_light/header['EXPTIME'] #converting to counts/s

    
    return target_light, bkg_light, noise_spectrum
