#https://juliaimagerecon.github.io/Sinograms.jl/stable/generated/examples/03-parallel-beam/
using Sinograms: SinoPar,  plan_fbp, fbp, fbp_sino_filter
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
using Sinograms: CtFanArc, CtFanFlat,CtFan # CtPar
using Sinograms:  plan_fbp, Window, Hamming, fdk, ct_geom_plot3
using ImageGeoms: ImageGeom, MaskCircle, fovs
using ImagePhantoms: ellipsoid_parameters, ellipsoid
using ImagePhantoms: radon, phantom
using MIRTjim: jim, prompt
using ImagePhantoms: Ellipsoid, ellipsoid
using ImagePhantoms: Cuboid, cuboid
using Sinograms: project_bdd, backproject_bdd
using Base.Threads: @threads
using Dates
using HDF5
# using RadonKA


"""
according to chat gpt 
ns::Int64 888: Number of detector elements along the s-axis (detector width).
nt::Int64 64: Number of detector elements along the t-axis (detector height).
ds::Float64 1.0239: Detector element size along the s-axis.
dt::Float64 1.0964: Detector element size along the t-axis.
offset_s::Float32 1.25: Offset of the detector elements along the s-axis.
offset_t::Float32 0.0: Offset of the detector elements along the t-axis.
na::Int64 984: Number of projection angles.
orbit::Float32 360.0: Total rotation angle of the source around the object (in degrees).
orbit_start::Float32 0.0: Starting angle of the source rotation.
source_offset::Float64 0.0: Offset of the source from the center of rotation.
dsd::Float64 949.075: Distance from the source to the detector.
dod::Float64 408.075: Distance from the object to the detector.
src::CtSourceCircle CtSourceCircle(): Represents the source trajectory, which in this case is a circular path.

Kowaluk Tomasz:
pierwsze 2 to po prostu rozdzielczość detektora, w tomografach jest na poziomie od 1000x1000 do 4000x4000
olejne 2 to rozmiar fizyczny piksela na detektorze i tu tez są różne wartości
offsety wynikają z kalibracji
liczba projekcji również zależy od systemu ja stosuję od 800 do 1500
zakres kątowy zawsze jest 360 a z dystansów źródło/detektor i obiekt/detektor liczymy powiększenie lub rozmiar woksela w zrekonstruowanym pliku.
W przypadku wiązki najczęściej jest to stożek a kąt zależy od źródła promieniowania. 
 



https://juliaimagerecon.github.io/Sinograms.jl/stable/generated/examples/07-fdk/
p = (ns = 130, ds = 0.3cm, nt = 80, dt = 0.4cm, na = 50, dsd = 200cm, dod = 40cm)
"""

# using Base.Threads: nthreads

# println("Number of threads available: ", nthreads())

# function threaded_map(f, collection)
#     results = Vector{eltype(f(collection[1]))}(undef, length(collection))
#     @threads for i in 1:length(collection)
#         results[i] = f(collection[i])
#     end
#     return results
# end

#was on more distant 147

"""
configuring computer tomography configuration
"""
function get_CTFAN_proj(ob,is_high_res=false)
    p = (ns = 300, ds = 0.15cm, nt = 300, dt = 0.15cm, na = 300, dsd = 250cm, dod = 100cm)
    # p = (ns = 800, ds = 0.15cm, nt = 800, dt = 0.15cm, na = 800, dsd = 200cm, dod = 40cm)
    # p = (ns = 800, ds = 0.15cm, nt = 800, dt = 0.2cm, na = 800, dsd = 200cm, dod = 40cm)
    # p = (ns = 500, ds = 0.2cm, nt = 500, dt = 0.2cm, na = 300, dsd = 40cm, dod = 40cm)
    if(is_high_res)
        p = (ns = 650, ds = 0.15cm, nt = 650, dt = 0.15cm, na = 650, dsd = 250cm, dod = 100cm)
    end
    
    rg = CtFanArc( ; p...)
    rayss=rays(rg)
    curr_time = Dates.now()
    proj_arc=zeros(Float32, p.ns*p.nt*p.na) #initialization
    radon(rayss, ob,proj_arc)
    # proj_arc=RadonKA.radon(collect(rayss), ob)
    proj_arc=reshape(proj_arc, (p.ns,p.nt, p.na))
    println("Time to get radon: ", Dates.now() - curr_time)
    return proj_arc,rg
end    

# function get_CTFAN_proj(ob,f_temp)
#     # p = (ns = 130, ds = 0.3cm, nt = 133, dt = 0.4cm, na = 207, dsd = 200cm, dod = 40cm)
#     p = (ns = 1000, ds = 0.3cm, nt = 1000, dt = 0.3cm, na = 984, dsd = 20cm, dod = 40cm)
#     dset = d_create(f_temp, "radon", Float32, (length,))
#     rg = CtFanArc( ; p...)
#     rayss=rays(rg)
#     curr_time = Dates.now()
#     proj_arc = radon(rayss, ob,f)
#     # proj_arc=RadonKA.radon(collect(rayss), ob)
#     proj_arc=reshape(proj_arc, (p.ns,p.nt, p.na))
#     println("Time to get radon: ", Dates.now() - curr_time)
#     return proj_arc,rg
# end    


function get_attenuation_coeficient()
    μ = 0.1 / cm # typical linear attenuation coefficient
    return μ
end
