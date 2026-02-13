#https://juliaimagerecon.github.io/Sinograms.jl/stable/generated/examples/03-parallel-beam/
using Sinograms: SinoPar, rays, plan_fbp, fbp, fbp_sino_filter
using ImageGeoms: ImageGeom, fovs, MaskCircle
using ImagePhantoms: SheppLogan, shepp_logan, radon, phantom
using Unitful: mm
using MIRTjim: jim, prompt
using ImagePhantoms: Object, phantom, radon, spectrum
using ImagePhantoms: Cylinder, cylinder
import ImagePhantoms as IP
using ImageGeoms: ImageGeom, axesf
using MIRTjim: jim, prompt, mid3
using FFTW: fft, fftshift, ifftshift
using LazyGrids: ndgrid
using Unitful: mm, unit, °
using Plots: plot, plot!, scatter!, default
using Plots # gif @animate
using PyCall
import Unitful
using Unitful
using Plots: plot, gui
using Unitful: cm
using Sinograms: CtFanArc, CtFanFlat # CtPar
using Sinograms: rays, plan_fbp, Window, Hamming, fdk, ct_geom_plot3
using ImageGeoms: ImageGeom, MaskCircle, fovs
using ImagePhantoms: ellipsoid_parameters, ellipsoid,half_sphere_x,half_sphere_y,half_sphere_z,half_sphere_z_b
using ImagePhantoms: radon, phantom
using MIRTjim: jim, prompt
using ImagePhantoms: Ellipsoid, ellipsoid
using ImagePhantoms: Cuboid, cuboid,cylinder_irr
import ImagePhantoms as IP




"""
construction is mainly based on /workspaces/synthethic_tomo/data/obraz.png
looking from the bottom we have a stem from aluminium with insulator poliethylene and copper 
then we get main part built from graphite that has multiple things in it
in the center it has central electrode which at the bottom is build from center - copper then poliethylene then insulator
asi it is in stem but at the top it is just copper (copper rod is longer then poliethylene and gets into 
top part of central elecrode that is build from graphite) around graphite part of central electrode 
there is air

From geometric perspective all is build from cylinders  
1) widest - is the main part of the chamber build from graphite together with a stem it creates all outer layer - its top and sides should be left not 
    covered by any other cylinder as it represents walls - through the bottom elements of electrode go through to the stem 
2) stem - is the part that is build from aluminium - it is the bottom part of the chamber it should have the same diameter as the widest 
    but it should be shorter (between 0.1 and 0.3 times the length of th widest) then the widest and it has aluminiuum density associated 
3) electrode - it is also cylinder and has multiple parts rom graphite, copper, poliethylene and insulator  
    External outline has 2 parts a bit bottom that is wider (from 20 to 10 % wider) and shorter (from 10 to 30 % shorter) then the top thinner part 
    top thinner part is entirely graphite and is surrounded by air (in order to get density of air one need to subtract the values from the 
    widest and then add them back in this top part) elongation of the top part gets to the top 10-15 % of bottom part is also build from graphite but 
    has the center with copper wire 
    copper wire has highest density of all and is in the center of the long axis of the widest and stem and electrode (they have common long axis) 
    is gets through all stem and all wider part of electrode  
    poliethylene is around the wider part of electrone on all of its sides and has lower density then copper under the poliethylene there is thin strip ofa luminium also 
    around whole center of the thicker part of the electrode  
    copper wire is surrounded with relatively thick insulator from poliethylene that is around 1.5 times wide then copper wire but do not reach the top of the wire 
    as the top of the wire is in the graphite part of the electrode 
Technical details 
 1)All of the cylinders has common long axis 
 2)the bottom wider part of electrode needs to be divided into stem part (b) and graphite part - as in practice we need to add or subtract 
 the density of the part we are on to get desired density - hence when we overlay copper wire on aluminium stem we need to add less then in case of overlying 
 copper wire on graphite part of the electrode
    
1) Current issues no aluminum strip below insulator is added
2) air and top thinner part of electrode is mis aligned in z axis
3) lower part of electrode is misaligned in z axis
4) parts from the same material has inconsistent density

"""
function ionic_chamber_p(
    center_main_cylinder
    ,main_cylinder_size
    ,graphite_top_thickness
    ,air_top_thickness
    ,air_thicknes
    ,copper_el_size
    ,polyethylene_size_dif
    ,insulator_thicknes
    ,central_electrode_size
    ,stem_height
    ,angle
    )

    polyethylene_size=(copper_el_size[1]+polyethylene_size_dif[1],copper_el_size[2]+polyethylene_size_dif[2],copper_el_size[3]+polyethylene_size_dif[3])
    #added in oder to mask negative half sphere
    widest=cylinder( (center_main_cylinder[1])cm, (center_main_cylinder[2])cm, (center_main_cylinder[3])cm
    ,(main_cylinder_size[1])cm, (main_cylinder_size[2])cm, (main_cylinder_size[3])cm,angle, 0., 0, 1.0f0)

    air=cylinder( (center_main_cylinder[1])cm, (center_main_cylinder[2])cm, (center_main_cylinder[3]+graphite_top_thickness)cm
    ,(main_cylinder_size[1]-air_thicknes)cm, (main_cylinder_size[2]-air_thicknes)cm, (central_electrode_size[3]+air_top_thickness)cm,angle, 0., 0, -1.0f0)

    # center_z_central=(center_main_cylinder[3]-(main_cylinder_size[3]/2)+graphite_top_thickness-air_top_thickness-central_electrode_size[3]*0.75)
    center_z_central=(center_main_cylinder[3]+central_electrode_size[3]*0.5)+air_top_thickness+graphite_top_thickness

    central_el=cylinder( (center_main_cylinder[1])cm, (center_main_cylinder[2])cm, (center_z_central)cm
    ,(central_electrode_size[1])cm, (central_electrode_size[2])cm, (central_electrode_size[3])cm,angle, 0., 0, 1.0f0)    
    # central_el=cylinder( (center_main_cylinder[1])cm, (center_main_cylinder[2])cm, (center_main_cylinder[3]+(main_cylinder_size[3]/2)-graphite_top_thickness-air_top_thickness-central_electrode_size[3]/2)cm
    # ,(central_electrode_size[1])cm, (central_electrode_size[2])cm, (central_electrode_size[3])cm,angle, 0., 0, 1.0f0)    

    copper_el=cylinder( (center_main_cylinder[1])cm, (center_main_cylinder[2])cm, (center_main_cylinder[3]+((main_cylinder_size[3]/2)+(copper_el_size[3]/2))-main_cylinder_size[3])cm
    ,(copper_el_size[1])cm, (copper_el_size[2])cm, (copper_el_size[3])cm,angle, 0., 0, 8.0f0)

    poliethylene=cylinder( (center_main_cylinder[1])cm, (center_main_cylinder[2])cm, (center_main_cylinder[3]+((main_cylinder_size[3]/2)+(polyethylene_size[3]/2))-main_cylinder_size[3])cm
    ,(polyethylene_size[1])cm, (polyethylene_size[2])cm, (polyethylene_size[3])cm,angle, 0., 0, 3.0f0)

    insulator=cylinder( (center_main_cylinder[1])cm, (center_main_cylinder[2])cm, (center_main_cylinder[3]+((main_cylinder_size[3]/2)+(copper_el_size[3]/2))-main_cylinder_size[3])cm
    ,(polyethylene_size[1]+insulator_thicknes)cm, (polyethylene_size[2]+insulator_thicknes)cm, (copper_el_size[3])cm,angle, 0., 0, -0.5f0)
    
    copper_el_b=cylinder( (center_main_cylinder[1])cm, (center_main_cylinder[2])cm, (center_main_cylinder[3]+((main_cylinder_size[3]/2)+(polyethylene_size[3]/2))-main_cylinder_size[3]-(stem_height*0.9))cm
    ,(copper_el_size[1])cm, (copper_el_size[2])cm, (stem_height)cm,angle, 0., 0, 7.5f0)

    poliethylene_b=cylinder( (center_main_cylinder[1])cm, (center_main_cylinder[2])cm, (center_main_cylinder[3]+((main_cylinder_size[3]/2)+(polyethylene_size[3]/2))-main_cylinder_size[3]-(stem_height*0.9))cm
    ,(polyethylene_size[1])cm, (polyethylene_size[2])cm, (stem_height)cm,angle, 0., 0, 2.5f0)

    insulator_b=cylinder( (center_main_cylinder[1])cm, (center_main_cylinder[2])cm, (center_main_cylinder[3]+((main_cylinder_size[3]/2)+(polyethylene_size[3]/2))-main_cylinder_size[3]-(stem_height*0.9))cm
    ,(polyethylene_size[1]+insulator_thicknes)cm, (polyethylene_size[2]+insulator_thicknes)cm, (stem_height)cm,angle, 0., 0, 0.0f0)
    
    aluminium_stem=cylinder( (center_main_cylinder[1])cm, (center_main_cylinder[2])cm, (center_main_cylinder[3]-(main_cylinder_size[3]/2)-(stem_height/2) )cm
    ,(main_cylinder_size[1])cm, (main_cylinder_size[2])cm, (stem_height)cm,angle, 0., 0, 0.5f0)    

    ob=[widest,air,central_el,copper_el,poliethylene,insulator,aluminium_stem,copper_el_b,poliethylene_b,insulator_b]
    
    vols=Dict("graphite_outer" => IP.volume(widest) - IP.volume(air)-IP.volume(poliethylene) ,
        "central_elecrode" => IP.volume(central_el),
        "aluminium_stem" => IP.volume(aluminium_stem)-IP.volume(insulator_b),
        "copper_electrode" => IP.volume(copper_el)+IP.volume(copper_el_b),
        "insulator" => (IP.volume(insulator)+IP.volume(insulator_b))-(IP.volume(poliethylene)+IP.volume(poliethylene_b)),
        "poliethylene" => (IP.volume(poliethylene)+IP.volume(poliethylene_b))-(IP.volume(copper_el)+IP.volume(copper_el_b))    
    )

    return ob,vols
end  

function ionic_chamber()
    center_main_cylinder = (0, 0, 0)
    main_cylinder_size=(4,4,9)
    graphite_top_thickness=1.5
    air_top_thickness=1.0
    air_thicknes=1.0
    copper_el_size=(0.4,0.4,2)
    polyethylene_size_dif=(1.1,1.1,-0.5)
    insulator_thicknes=0.5
    central_electrode_size=(1,1,3)
    stem_height=2.0
    angle=0

    return ionic_chamber_p(center_main_cylinder,main_cylinder_size,graphite_top_thickness,air_top_thickness,air_thicknes,copper_el_size,polyethylene_size_dif
    ,insulator_thicknes,central_electrode_size,stem_height,angle)
        
end  


function empty_cylinder_with_half_sphere_bottom()
    center_cylinder = (0, 0, 0)
    bigger_cyl_size=(2,4,8)
    cylinder_wall_thicness=0.5
    cylinder_bottom_curvature=1.2
    cylinder_top_curvature=0.9
    angle=π/6
    density_inside=0.15

    
    return empty_cylinder_with_half_sphere_bottom_p(center_cylinder,bigger_cyl_size,cylinder_wall_thicness
            ,cylinder_bottom_curvature,cylinder_top_curvature,angle,density_inside
            ,1.0,(0.4,0.4),2.5
            ,1.0,(0.4,0.4),1.5)

end  

function volume_of_elliptical_cylinder(pipe_cross_section, overlay_length)
    a = pipe_cross_section[1] / 2
    b = pipe_cross_section[2] / 2
    h = overlay_length
    volume = π * a * b * h
    return volume
end

function empty_cylinder_with_half_sphere_bottom_p(
     center_cylinder #needs to be within the max(bigger cylinder size) from edges of the image
    ,bigger_cyl_size # its max need to be less then 0.75 of image size
    ,cylinder_wall_thicness# less then 1/4 of the bigger_cyl_size more then 0.1
    ,cylinder_bottom_curvature# min 0 max half of the cylinder size[3]
    ,cylinder_top_curvature# min 0 max half of the cylinder size[3]
    ,angle #between 0 and 2 pi
    ,density_inside#between0.05 and 0.9 
    ,pipe_len #0.2 to 0.8 of the cylinder size[3]
    ,pipe_cross_section #0.02 to 0.1 of the cylinder size[1] and cylinder size[2]
    ,pipe_density #0.1 to 0.9
    ,dispenser_len #0.05 to 0.2 of the cylinder size[3]
    ,dispenser_cross_section #0.02 to 0.1 of the cylinder size[1] and cylinder size[2]
    ,dispenser_density #0.1 to 0.5
    ,len_cut
    ,menisc_radius
    ,dual_phase_percentage
    ,density_inside_b
     )

    dispenser_cross_section=(pipe_cross_section[1]+dispenser_cross_section[1],pipe_cross_section[2]+dispenser_cross_section[2])
    halph_s_bigger_size=(bigger_cyl_size[1]-cylinder_wall_thicness,bigger_cyl_size[2]-cylinder_wall_thicness,cylinder_bottom_curvature)
    halph_s_bigger_size_top=(bigger_cyl_size[1]-cylinder_wall_thicness,bigger_cyl_size[2]-cylinder_wall_thicness,cylinder_top_curvature)

    # rrr=((halph_s_bigger_size_top[1]-cylinder_wall_thicness)cm, (halph_s_bigger_size_top[2]-cylinder_wall_thicness)cm, (halph_s_bigger_size_top[3])cm)
    # print("\n aaaaa cylinder_wall_thicness $cylinder_wall_thicness bigger_cyl_size $bigger_cyl_size halph_s_bigger_size_top $halph_s_bigger_size_top diff $rrr  halph_s_bigger_size $halph_s_bigger_size \n")

    if(dual_phase_percentage==1.0)
        ob2_a=cylinder( center_cylinder[1]cm, center_cylinder[2]cm,(center_cylinder[3])cm, (bigger_cyl_size[1]-(cylinder_wall_thicness))cm, (bigger_cyl_size[2]-cylinder_wall_thicness)cm, (bigger_cyl_size[3])cm, angle, 0, 0, (-(1-density_inside))/2)
        ob2_b=cylinder( center_cylinder[1]cm, center_cylinder[2]cm,(center_cylinder[3])cm, (bigger_cyl_size[1]-(cylinder_wall_thicness))cm, (bigger_cyl_size[2]-cylinder_wall_thicness)cm, (bigger_cyl_size[3])cm, angle, 0, 0, (-(1-density_inside))/2)
    else
        len1=bigger_cyl_size[3]*dual_phase_percentage
        len2=bigger_cyl_size[3]-len1
        L=len1+len2
        z_center=center_cylinder[3]
        z1 = z_center - (L / 2) + (len1 / 2)
        z2 = z_center + (L / 2) - (len2 / 2)

        ob2_a=cylinder( center_cylinder[1]cm, center_cylinder[2]cm,(z1)cm
        , (bigger_cyl_size[1]-(cylinder_wall_thicness))cm
        , (bigger_cyl_size[2]-cylinder_wall_thicness)cm
        , (len1)cm, angle, 0, 0, -(1-density_inside))

        ob2_b=cylinder( center_cylinder[1]cm, center_cylinder[2]cm,(z2)cm
        , (bigger_cyl_size[1]-(cylinder_wall_thicness))cm
        , (bigger_cyl_size[2]-cylinder_wall_thicness)cm
        , (len2)cm, angle, 0, 0, -(1-density_inside_b))

    end    

        
    ob_cut=cylinder( center_cylinder[1]cm, center_cylinder[2]cm,((center_cylinder[3]+((bigger_cyl_size[3]/2)- (len_cut/2))))cm, ((bigger_cyl_size[1]-(cylinder_wall_thicness))*0.99)cm
    , ((bigger_cyl_size[2]-cylinder_wall_thicness)*0.99)cm, (len_cut)cm, angle, 0, 0,(-density_inside))
    ob_menisc_cut=half_sphere_z_b( center_cylinder[1]cm, center_cylinder[2]cm ,((center_cylinder[3]+((bigger_cyl_size[3]/2))))cm
    , (bigger_cyl_size[1]*1.5)cm, (bigger_cyl_size[2]*(1.5))cm, (menisc_radius)cm, angle, 0, 0,1.0)

    ob3=cylinder( center_cylinder[1]cm, center_cylinder[2]cm,center_cylinder[3]cm, bigger_cyl_size[1]cm, bigger_cyl_size[2]cm, bigger_cyl_size[3]cm,angle, 0, 0, 1.0f0)
    ob_cyl_mask=cylinder( center_cylinder[1]cm, center_cylinder[2]cm,center_cylinder[3]cm, (bigger_cyl_size[1]*1.05)cm, (bigger_cyl_size[2]*1.05)cm, bigger_cyl_size[3]cm,angle, 0, 0, 1.0f0)


    ob4=half_sphere_z((center_cylinder[1])cm, (center_cylinder[2])cm,
     ((center_cylinder[3]-bigger_cyl_size[3]/2))cm, #z
     (halph_s_bigger_size[1]-cylinder_wall_thicness)cm, (halph_s_bigger_size[2]-cylinder_wall_thicness)cm, (halph_s_bigger_size[3]-(cylinder_wall_thicness))cm, angle, 0., 0, -2+density_inside)
    bottom_cyl_start=((center_cylinder[3]-bigger_cyl_size[3]/2))+cylinder_wall_thicness
    ob5=half_sphere_z((center_cylinder[1])cm, (center_cylinder[2])cm
    , (((center_cylinder[3]-bigger_cyl_size[3]/2)))cm, halph_s_bigger_size[1]cm, halph_s_bigger_size[2]cm, halph_s_bigger_size[3]cm, angle, 0., 0, 1.0f0)
     
    ob4b=half_sphere_z((center_cylinder[1])cm, (center_cylinder[2])cm, ((center_cylinder[3]+bigger_cyl_size[3]/2))cm, (halph_s_bigger_size_top[1])cm, (halph_s_bigger_size_top[2])cm, (halph_s_bigger_size_top[3])cm, angle, 0., 0, -(1))
    ob5b=half_sphere_z((center_cylinder[1])cm, (center_cylinder[2])cm, ((center_cylinder[3]+bigger_cyl_size[3]/2))cm, (halph_s_bigger_size_top[1]+cylinder_wall_thicness)cm, (halph_s_bigger_size_top[2]+cylinder_wall_thicness)cm, (halph_s_bigger_size_top[3]+cylinder_wall_thicness)cm, angle, 0., 0, 1.0f0)
    
    ob5b_mask=half_sphere_z((center_cylinder[1])cm, (center_cylinder[2])cm, ((center_cylinder[3]+bigger_cyl_size[3]/2))cm, 
    ((halph_s_bigger_size_top[1]+cylinder_wall_thicness)*1.05)cm, ((halph_s_bigger_size_top[2]+cylinder_wall_thicness)*1.05)cm, ((halph_s_bigger_size_top[3]+cylinder_wall_thicness)*1.05)cm, angle, 0., 0, 1.0f0)
    
    #added in order to mask negative half sphere
    # ob6=cylinder( (center_cylinder[1])cm, (center_cylinder[2])cm, ((bottom_cyl_start-(cylinder_wall_thicness)  ))cm, (bigger_cyl_size[1]-cylinder_wall_thicness)cm, (bigger_cyl_size[2]-cylinder_wall_thicness)cm, (cylinder_wall_thicness)cm,angle, 0., 0, 1.0f0)
    # ob6=cylinder( (center_cylinder[1])cm, (center_cylinder[2])cm, ((center_cylinder[3]-bigger_cyl_size[3]/2))cm, (halph_s_bigger_size[1]-cylinder_wall_thicness)cm, (halph_s_bigger_size[2]-cylinder_wall_thicness)cm, (cylinder_wall_thicness)cm,angle, 0., 0, 1.0f0)
    pipe_start=center_cylinder[3]+(bigger_cyl_size[3]/2)+cylinder_top_curvature-(pipe_len/2)
    pipe=cylinder( center_cylinder[1]cm, center_cylinder[2]cm,(pipe_start)cm, pipe_cross_section[1]cm, pipe_cross_section[2]cm, (pipe_len)cm,angle, 0, 0, pipe_density)
    pipe_in=cylinder( center_cylinder[1]cm, center_cylinder[2]cm,(pipe_start)cm, pipe_cross_section[1]cm/2, pipe_cross_section[2]cm/2, (pipe_len)cm,angle, 0, 0, (-pipe_density))
    dispenser_start=center_cylinder[3]+(bigger_cyl_size[3]/2)+cylinder_top_curvature#-(dispenser_len)
    plastic_dispenser=cylinder( center_cylinder[1]cm, center_cylinder[2]cm,(dispenser_start)cm, dispenser_cross_section[1]cm, dispenser_cross_section[2]cm, (dispenser_len)cm,angle, 0, 0, dispenser_density)

    # ob=[ob2,ob3,ob4,ob5,ob4b,ob5b,pipe,pipe_in,plastic_dispenser]#ob6
    ob=[ob2_a,ob2_b,ob3,pipe,pipe_in,plastic_dispenser,ob_cut]#ob6

    vol1=IP.volume(ob2_a)+IP.volume(ob2_b)
    vol2=IP.volume(ob5)
    # vol3=IP.volume(ob5b)
    vol_pipe=IP.volume(pipe)
    #we need to calculate what percentege of the pipe is in the "fluid" so that is overlaying the ob2
    # overlay is important for us just in z dimension  
    pipe_center_z=pipe_start
    pipe_beg_z=pipe_start-(pipe_len/2)
    pipe_end_z=pipe_start+(pipe_len/2)

    ob2_beg_z=center_cylinder[3]-(bigger_cyl_size[3]/2)
    ob2_end_z=center_cylinder[3]+(bigger_cyl_size[3]/2)

    overlay=0
    if pipe_beg_z<=ob2_end_z && pipe_end_z>=ob2_beg_z
        overlay=abs(min(pipe_end_z,ob2_end_z)-max(pipe_beg_z,ob2_beg_z))
        # overlay=abs(max(pipe_end_z,ob2_end_z)-min(pipe_beg_z,ob2_beg_z))
        # overlay=abs(pipe_end_z-ob2_end_z)
    end
    
    #get a volume of a cylinder of radius of the elipsoid that is crossection (pipe_cross_section[1]/2) and (pipe_cross_section[2]/2) and a length of overlay
    whole_pipe=cylinder( center_cylinder[1]cm, center_cylinder[2]cm,
    (pipe_center_z-overlay/2)cm
    , pipe_cross_section[1]cm, pipe_cross_section[2]cm, (overlay)cm,angle, 0, 0, pipe_density)
    whole_pipe_volume=IP.volume(whole_pipe)

    # whole_pipe_volume=volume_of_elliptical_cylinder(pipe_cross_section, overlay)
    # inside_of_pipe_volume= volume_of_elliptical_cylinder((pipe_cross_section[1]/2,pipe_cross_section[2]/2), overlay)
    inside_of_pipe_volume= IP.volume(cylinder( center_cylinder[1]cm, center_cylinder[2]cm,(pipe_end_z)cm, pipe_cross_section[1]cm/2, pipe_cross_section[2]cm/2, (overlay)cm,angle, 0, 0, (-pipe_density)))
    # print("\n for debug should equal $(volume_of_elliptical_cylinder(pipe_cross_section, pipe_len)) b: $(vol_pipe)")
    res_vol=vol1-vol2-(whole_pipe_volume)+inside_of_pipe_volume-IP.volume(ob_cut)

    return ob,(Dict("can_inside"=>(res_vol) )),ob4,ob5,ob4b,ob5b,ob_cut,ob_menisc_cut,ob_cyl_mask,ob5b_mask

end  


function get_geometric_onjects()
    # return ionic_chamber()
    return empty_cylinder_with_half_sphere_bottom()

end

#TODO get segmentation - render separately each part that we are intrested without noise and then fuse all voxels with the same value
