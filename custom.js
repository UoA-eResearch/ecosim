/**
 * Returns a number whose value is limited to the given range.
 *
 * Example: limit the output of this computation to between 0 and 255
 * (x * 255).clamp(0, 255)
 *
 * @param {Number} min The lower boundary of the output range
 * @param {Number} max The upper boundary of the output range
 * @returns A number in the range [min, max]
 * @type Number
 */
Number.prototype.clamp = function(min, max) {
    return Math.min(Math.max(this, min), max);
};

function HN() {
    B=400;           //biomass mgC/m2
    B_mac=300;     //biomass mgC/m2
    D_mac=10.0;     //detritus 
    D_MPB=100.0;     //detritus 
    N=100.0;     //nutrient
}

function LN() {
    B=1;           //biomass mgC/m2
    B_mac=1;     //biomass mgC/m2
    D_mac=0;     //detritus 
    D_MPB=0;     //detritus 
    N=0;     //nutrient
}

function init() {

    depth=1; //what it says, meters

    mortB=0.01;         //mortality of mpb per timestep. Should be dynamic. 1/day (M&F, JGR 2012)
    growthrateB=10*0.0062;     //1/day, max growth of mpb, should depend at least on //mud and nutrients, units: 1/day from Niathany 
    grazingB=0.0002;       //1/day, ?????????? 
    K=1200;            //mgCha/m2carrying capacity, max mpb biomass. //Ecolog. Modelling 2016 Sinha 
    Xb=0;              //minimum background mpb biomass, mgCha/m2

    dfine=0.000012;	// diameter of fine sed, in m */
    dcoarse=0.000125;   // diameter of coarse sed, in m */
    rho=1000;             	// water density, kg/m^3 */
    rhoS=2600;           	// sediment density, kg/m^3 */
    //Shieldscriticalcoarse has been evaluated from the shields diagram
    //for dcoarse=0.000125, see also page 106 of Dynamics of marine sand, soulsby 1997*/
    Shieldscriticalcoarse=0.065;
    g=9.81;  //gravity m/s2

    growthB_mac_max=0.7;       //1/day, max growth of macomona. 
    K_mac=1200;            //mgC/m2carrying capacity, max number of bugs possible, should depend at least on d50 (simon's paper). 
    Xb_mac=0;              //minimum number of bugs, mgCha/m2
    wetW2Carbon=0.05*(10**-3)*10**6; //mg C per m2 per individual (from individuals to carbon)
    //assumig a shell lenth of 1cm, gives a volume of 0.01*0.01*0.01 m3
    //(10^-3 factor in eq below). We assume the specific weight of clams
    //is 1000kg/m^3 and 0.05 is the wet weight (AFDW in Ricciardi and
    //Bourget, 1998 meps) to carbon, convert from Kg to mg
    //**********CHECK DAN PRATT THESIS, PAGE 62

    mortB_mac=0.01;         //mortality of macomona per timestep. Should be dynamic. ind/day from judi spreadsheet ...
    //We actually realized that Judi's approach results in too much of a large
    //contribution and makes the population go negative so ...
    // we decided to have just a percentage of the macomona
    //mortB_mac=mortB_mac*wetW2Carbon; //mgC/m2ind

    //initial and boundary conditions
    // B=400*ones(tmax/100,1);           //biomass mgC/m2
    // B_mac=300*ones(tmax/100,1);     //biomass mgC/m2
    // D_mac=10.0*ones(tmax/100,1);     //detritus 
    // D_MPB=10.0*ones(tmax/100,1);     //detritus 
    // N=1000000.0*ones(tmax/100,1);     //nutrient 
    eval(scenario + "()")
    IN=1000;                   //external source of nutrient. It needs to have the same units as N & Kn
    density_mac=0.1; //this can only change between 0 and 1

    //Light effect on growth rate (it would be nice to link this to tidal range)
    I0=300;          //I0=300miuE/m2/s from M&F
    kd=0.4;    //changes depnding on turbidity (mariana) units: 1/meters
    I=I0*Math.exp(-kd*depth);
    Ki=20;          //Ki=20miuE/m2/s from Bouletreau et al. 2008..probably not for local species
    light=I/(I+Ki); //light reduction term

    //Mud content effect on growth rate (set to 1 for now. The growth should
    //actually depend on nutrient availability which is in turn a function of
    //the mud content)
    // mudpc=[0.1:0.1:1];
    //mudpcf=0.9;
    //mud_content=2*(1-1./(1+exp(-mudpcf))); //controls the flux of nutrients out
    //of the sediment... not sure how we came up with that
    mud_content=0.1;
    //This is a pore water efficiency factor. More flushing with sand, less
    //flushing if it's all mud. Based on hydraulic conductivity?

    Kn=10;      //micrograms of Nitrogen/l, from Hydrobiologia, Naithani 2016. Biological uptake, species-dependent
    coeff_BD=0.16;  //original value Naithani et al 2016, 1/day
    conversion_C2N=5.88; //Naithani et al.
    depthbenthiclayer=0.01; //Ecolog. Modelling 2016 Sinha 
    //K_uptakenutrient=1.6; //same unit as nutrents mmol*m^-2 day^-1. Data from Sundback & Granelli
    K_uptakenutrient=50; //same unit as nutrents mmol*m^-2 day^-1. Data from Sundback & Granelli
    coefflushing=0.00001;

    //initialization
    B_new=B;B_mac_new=B_mac;
    D_MPB_new=D_MPB; D_mac_new=D_mac; N_new=N;N_new_mac=N/2;N_new_MPB=N/2;
    dt=0.00001;
    dtsave=1000;
    dt = 1;
    updateControls();
    updateBoxes();

    var N_new_mac_trace = {
        y: [N_new_mac],
        customdata: "N_new_mac",
        mode: 'lines',
        name: '<i>Macomona liliana</i> nutrients mgC/m²',
        line: {
            color: "00ffffff"
        },
        hoverlabel: {
            namelength: -1
        },
        visible: "legendonly"
    };

    var N_new_MPB_trace = {
        y: [N_new_MPB],
        customdata: "N_new_MPB",
        mode: 'lines',
        name: 'Microphytobenthos nutrients mgC/m²',
        line: {
            color: "dd55ffff"
        },
        hoverlabel: {
            namelength: -1
        },
        visible: "legendonly"
    }

    var B_trace = {
        y: [B],
        customdata: "B",
        mode: 'lines',
        line: {
            color: "71c837ff"
        },
        hoverlabel: {
            namelength: -1
        },
        name: 'Microphytobenthos Biomass mgC/m²'
    };

    var B_mac_trace = {
        y: [B_mac],
        customdata: "B_mac",
        mode: 'lines',
        line: {
            color: "2a7fffff"
        },
        hoverlabel: {
            namelength: -1
        },
        name: '<i>Macomona liliana</i> Biomass mgC/m²'
    };

    var D_mac_trace = {
        y: [D_mac],
        customdata: "D_mac",
        mode: 'lines',
        line: {
            color: "ff9955ff"
        },
        hoverlabel: {
            namelength: -1
        },
        name: '<i>Macomona liliana</i> detritus mgC/m²',
        visible: "legendonly"
    };

    var D_MPB_trace = {
        y: [D_MPB],
        customdata: "D_MPB",
        mode: 'lines',
        name: 'Microphytobenthos detritus mgC/m²',
        line: {
            color: "ff5555ff"
        },
        hoverlabel: {
            namelength: -1
        },
        visible: "legendonly"
    };

    var mud_content_trace = {
        y: [mud_content],
        customdata: "mud_content",
        mode: 'lines',
        line: {
            color: "996633ff"
        },
        hoverlabel: {
            namelength: -1
        },
        name: 'Mud Content',
        visible: "legendonly"
    }

    var data = [ B_mac_trace, B_trace, N_new_mac_trace, N_new_MPB_trace, D_mac_trace, D_MPB_trace, mud_content_trace ];

    var layout = {
        xaxis: {
          title: 'Days'
        },
    };

    Plotly.newPlot('plot', data, layout, {responsive: true});
}

$("#plot").on('plotly_hover', function(event, data){
    if (playing) return;
    svg = $("svg", $("#boxes")[0].contentDocument);
    if (svg.length == 0) return;

    for (var i in data.points) {
        var p = data.points[i];
        var name = p.data.customdata;
        var val = p.y;
        console.log(name, val);
        switch(name) {
            case "B_mac":
                scaleBox("#B_mac", val / 300);
                scaleBox("#B_mac_D_mac", val / 300);
                break;
            case "B":
                scaleBox("#B", val / 400);
                scaleBox("#B_D_MPB", val / 400);
                scaleBox("#B_B_mac", val / 400);
                break;
            case "N_new_mac":
                scaleBox("#nutrients_depth", val / 50);
                scaleBox("#nutrients_depth_surface_nutrients", val / 50);
                break;
            case "N_new_MPB":
                scaleBox("#surface_nutrients", val / 50);
                scaleBox("#surface_nutrients_B", val / 50);
                break;
            case "D_mac":
                scaleBox("#D_mac", val / 10);
                scaleBox("#D_mac_nutrients_depth", val / 10);
                break;
            case "D_MPB":
                scaleBox("#D_MPB", val / 100);
                scaleBox("#D_MPB_surface_nutrients", val / 100);
                break;
        }
    }
});

$("#plot").on("plotly_unhover", function() {
    updateBoxes();
});

function advance() {
    N_old=N_new; 
    N_old_mac=N_new_mac;
    N_old_MPB=N_new_MPB;
    B_old=B_new;B_mac_old=B_mac_new;
    D_MPB_old=D_MPB_new; D_mac_old=D_mac_new;
        
    //biomass dynamics
    //the reduced uptake when sediment is muddier. Look at 14.70, p. 479
    //of the big modelling book
    if (IN < N_old_MPB) {
        growthB=growthrateB*light*(N_old_MPB)/(K_uptakenutrient+(N_old_MPB)); //we need to check the uptake value
    } else {
        growthB=growthrateB*light*(IN)/(K_uptakenutrient+(IN)); //we need to check the uptake value  
    }
    growthB_mac=growthB_mac_max*B_old/K;
    //with carrying capacity
        
    //MACOMONA DYNAMICS
    K_mac_new=(21.6+0.146*(mud_content*100)-0.005*((mud_content*100)^2))*10000/353.4; // sIMON ET AL 2003 meps. uNITS ARE INDIVIDUALS/QUADRAT (0.25^2)

    B_new=B_old+dt*(growthB*B_old*(1-B_old/K)-mortB*(B_old-Xb)-growthB_mac*B_mac_old*(1-B_mac_old/K_mac_new));

    //without carrying capacity
    //B(i)=B(i-1)+dt*(growthB*B(i-1)-mortB*(B(i-1)-Xb)-growthB_mac*B_mac(i-1)*(1-B_mac(i-1)/K_mac(i-1)));

    //we need to fix this bit
    // K_mac_carb(i)=wetW2Carbon*K_mac(i-1); //mg C per m2 
    if (N_old_mac>333333) mortB_mac=0.1; 
    B_mac_new=B_mac_old+dt*(growthB_mac*B_mac_old*(1-B_mac_old/K_mac_new)-mortB_mac*B_mac_old);
    //D(i)=coeff_BD*mortB*B(i-1)-lostD*D(i-1)-convertedD*D(i-1);
    //coeff_BD could be different for mac and mpb

    //coefflushing=c1*mud_content+c2*B_mac_old;
    
    D_MPB_new=D_MPB_old+dt*(mortB*B_old-coeff_BD*D_MPB_old);
    D_mac_new=D_mac_old+dt*(mortB_mac*B_mac_old-coeff_BD*D_mac_old);
    N_new_MPB=N_old_MPB+dt*(coeff_BD*D_MPB_old-growthB*B_old*(1-B_old/K))/depthbenthiclayer/conversion_C2N+dt*N_old_mac*coefflushing; 
    N_new_mac=N_old_mac+dt*(coeff_BD*D_mac_old)/depthbenthiclayer/conversion_C2N-dt*N_old_mac*coefflushing; 
    N_new=N_new_MPB+N_new_mac;
    
    K_mac=K_mac_new;
    B=B_new;
    B_mac=B_mac_new;
    D_MPB=D_MPB_new;
    D_mac=D_mac_new;
    N=N_new; 
    N_MPB=N_new_MPB;
    N_mac=N_new_mac;
    // [N(i-1) coeff_BD*D_MPB(i-1) coeff_BD*D_mac(i-1) -growthB*B(i-1)]
    // return
    if (N_new<0) {
        N_new = 0;
        console.warn("no more N");
    }
    if (N_new_MPB<0) {
      N_new_MPB=0;
      console.warn("no more N_new_MPB");
    }
    if (N_new_mac<0) {
      N_new_mac=0;
      console.warn("no more N_new_mac");
    }
    //****WARNING !!!!! WE NEED TO CHECK
    //depthbenthiclayer/conversion_C2N (28/8/2018)
    //MPB_n = B/conversion_C2N/depthbenthiclayer;
    //Detritus_mac_n = D_mac/conversion_C2N/depthbenthiclayer;
    //Detritus_MPB_n = D_MPB/conversion_C2N/depthbenthiclayer;
    // plot N, B_mac, B, N_MPB, N_mac - Biomass mgC/m2
    updateControls();
    updateBoxes();
    Plotly.extendTraces('plot', {
        //B_mac_trace, B_trace, N_new_mac_trace, N_new_MPB_trace, D_mac_trace, D_MPB_trace, mud_content_trace
        y: [[B_mac], [B], [N_new_mac], [N_new_MPB], [D_mac], [D_MPB], [mud_content]]
    }, [0,1,2,3,4,5,6])
}

function scaleBox(boxID, scale) {
    var box = $(boxID, svg);
    if (box.length == 0) {
        console.warn(boxID + " isn't ready yet!");
        return;
    }
    var bbox = box[0].getBBox();
    var cx = bbox.x + bbox.width / 2;
    var cy = bbox.y + bbox.height / 2;
    box.css({
        "transform-origin": cx + "px " + cy + "px",
        "transform": "scale(" + scale.clamp(.5, 1.2) + ")"
    });
}

function updateBoxes() {
    svg = $("svg", $("#boxes")[0].contentDocument);
    if (svg.length == 0) return;

    scaleBox("#B_mac", B_mac_new / 300);
    scaleBox("#B_mac_D_mac", B_mac_new / 300);
    scaleBox("#B", B_new / 400);
    scaleBox("#B_D_MPB", B_new / 400);
    scaleBox("#B_B_mac", B_new / 400);
    scaleBox("#nutrients_depth", N_new_mac / 50);
    scaleBox("#nutrients_depth_surface_nutrients", N_new_mac / 50);
    scaleBox("#surface_nutrients", N_new_MPB / 50);
    scaleBox("#surface_nutrients_B", N_new_MPB / 50);
    scaleBox("#D_mac", D_mac_new / 10);
    scaleBox("#D_mac_nutrients_depth", D_mac_new / 10);
    scaleBox("#D_MPB", D_MPB_new / 100);
    scaleBox("#D_MPB_surface_nutrients", D_MPB_new / 100);
}

function updateControls() {
    $("input").each(function() {
        var val = eval(this.id);
        if (this.id != "mud_content") {
            val = Math.round(val);
        }
        $(this).val(val);
    })
}

scenario = "HN";
init();

$("#reset").click(function () {
    init();
});

playing = false;

setInterval(function() {
    if (playing) {
        advance();
    }
}, 10);

$("#play").click(function() {
    playing = !playing;
    if (playing) {
        $("#play i").attr("class", "fas fa-pause");
        $("#play span").text("Pause");
        $("#step").hide();
    } else {
        $("#play i").attr("class", "fas fa-play");
        $("#play span").text("Play");
        $("#step").show();
    }
});

$("#step").click(function () {
    advance();
});

$("input").on("input", function(e) {
    console.log(this.id, this.value);
    eval(this.id + "=" + this.value);
    updateBoxes();
});

$("#scenario").change(function() {
    scenario = this.value;
    init();
});