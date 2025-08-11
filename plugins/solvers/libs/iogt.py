""" 
Default parameters for model setup and payoff metrics
    1. Changes here is cause issues with loading previous simulations, compatibility not assured
    2. These values may be modified by the solver as needed to obtain valid solutions
"""
#initializations
import os
import numpy as np
if __package__ is None or __package__ == '':
    import properties
    import stress
    from units import *
else:
    from . import properties
    from . import stress
    from .units import *
water = properties.water()
deg=deg
yr=yr
cP=cP
darcy=darcy
mD=mD
g=g

#model definition
class setup:
    def __init__(self):
        #GeoDT solver parameters
        self.Strategy_Iterations_units = 10
        self.Strategy_DiscreteFractures_path = os.getcwd()+'\\dfn.csv'
        self.Strategy_Solver_type = 'gg' #'Gringarten' 'GeoDT' 'gg' 'gt'
        self.Strategy_Design_type = 'EGS' #'ALL' 'EGS' 'AGS' 'CGS' 'O&G' 'DCM' 'FGS'
        self.Strategy_Target_type = 'Deep' #'Deep','Temp'
        self.Strategy_Stim_type = 'Radi' #'Volu','Radi'
        self.Strategy_SecondOrderNoise_ratio = 0.0 #0.0 = no second order randomness, 1.0 = full, 2.0 = amplified #TODO
        self.Well_TargetTemp_K = 273.15+325 #K
        self.Well_TargetDepth_m = 6000.0 #m
        self.Well_ProducerCount_wells = 4 #wells
        self.Well_Pattern_count = 1 #pattern #TODO: needs to be implemented
        self.Well_Spacing_m = 400.0 #m
        self.Well_DeviatedLength_m = 800.0 #m
        self.Well_ProducerProportion_ratio = 0.6 #m/m
        self.Well_Azimuth_rad = 90.0*deg #rad
        self.Well_Dip_rad = 30.0*deg #rad
        self.Well_RotationPhase_rad = 0.0*deg #rad
        self.Well_RotationToe_rad = 0.0*deg #rad
        self.Well_RotationSkew_rad = 0.0*deg #rad
        self.Well_ProductionDiameter_m = 0.0254*8.625 #m #TODO: Implement complete well plan
        self.Well_WallThickness_m = 0.0254*0.75 #m #TODO: Implement complete well plan
        self.Well_AnnularThickness_m = 0.0254*1.25 #m #TODO: Implement complete well plan
        self.Well_SegmentCount_units = 10.0 #units, number of breaks in the well for AGS, first is conductor #TODO: Implement complete well plan
        self.Well_Roughness_metric = 80.0
        self.Circulation_InjectionRate_m3ps = 0.10 #m3/s
        self.Stimulation_InjectionRate_m3ps = 0.4 #m3/s
        self.Stimulation_TargetVolume_m3 = 50000.0 #m3
        self.Stimulation_TargetRadius_m = 525.0 #m3
        self.Domain_Size_m = 1000.0 #m
        self.Rock_ThermalGradient_Kpm = 0.060 #K/m
        self.Rock_Density_kgpm3 = 2700.0 # kg/m3
        self.Rock_ThermalConductivity_WpmK = 2.5 # W/m-K
        self.Rock_HeatCapacity_kJpm3K = 2063.0 # kJ/m3-K
        self.Rock_Biot_ratio = 0.0 #Pa/Pa
        self.Air_Temperature_K = 273.15 # K
        self.Air_Pressure_Pa = 0.101*1e6 #Pa
        self.Rock_YoungsModulus_Pa = 50.0*1e9 #Pa
        self.Rock_PoissonRatio_ratio = 0.3 #
        self.Stress_ShminToSv_ratio = 0.5 #Pa/Pa
        self.Stress_ShmaxToShmin_ratio = 1.5 #Pa/Pa
        self.Stress_ShminAzimuth_rad = 110.0*deg #rad
        self.Stress_ShminDip_rad = 30.0*deg #rad
        self.Stress_Variance_rad = 0.5*deg #rad
        self.Shear_SlipLength_ratio = 10.0**-2.0 #m/m #TODO
        self.Shear_DilationSlip_ratio = 0.200 #m/m #TODO
        self.Shear_HydraulicDilation_ratio = 0.6 #m/m #TODO
        self.Tension_FractureCompressibility_1pPa = 2.9e-8 #1/Pa #TODO
        self.Stimulation_ProppantCompressibility_1pPa = 2.9e-8 #1/Pa #TODO
        self.JointSets_InitialDilation_m = 0.00010 #m #TODO
        self.Domain_BoundaryAperture_m = 0.003 #m
        self.Fracture_Tortuosity_ratio = 0.80 #m/m #TODO
        self.Cement_ThermalConductivity_WpmK = 2.0 # W/m-K
        self.Cement_HeatCapacity_kJpm3K = 2000.0 # kJ/m3-K
        self.Power_GeneralEfficiency_ratio = 0.85 # kWe/kWt
        self.Power_LifeSpan_s = 30*yr #years
        self.Domain_TimeSteps_steps = 30 #steps
        self.Circulation_Backpressure_Pa = 1.0*1e6 #Pa
        self.Circulation_InjectionTemperature_K = 95.0+273.15 #K
        self.Power_ConvectionCoefficient_kWpm2K = 3.0 #kW/m2-K
        self.Rock_Permeability_m2 = 0.1*mD #m2
        self.Stimulation_ProppantPermeability_m2 = 100.0*darcy #m2 #TODO
        self.JointSets_Set1_fractures = 16 #fractures
        self.JointSets_Set2_fractures = 8 #fractures
        self.JointSets_Set3_fractures = 4 #fractures
        self.JointSets_Variance_rad = 7.0*deg #rad
        self.JointSets_DiameterMin_m = 10.0 #m
        self.JointSets_DiameterMax_m = 1000.0 #m
        self.Stimulation_Clusters_clusters = 3 #clusters
        self.Stimulation_PerfDiameter_m = 0.013 #m
        self.Stimulation_PerfPerCluster_perfs = 6 #holes
        self.Stimulation_ProppantConcentration_m3pm3 = 0.045 #m3/m3
        self.Domain_PressureIncrement_Pa = 0.1*1e6 #Pa
        self.Circulation_AllowableOverpressure_Pa = 0.99 #pinj/s3 #Pa
        self.Shear_FrictionCoefficient_ratio = 0.6 #T:N #TODO
        self.Shear_Cohesion_Pa = 1e6 #Pa #TODO
        self.Tension_Cohesion_Pa = 0.1*1e6 #Pa
        self.Tension_FrictionAngle_rad = 30.0*deg #rad
        self.Tension_Toughness_Pasqrtm = 1.5*1e6 #Pa-m**0.5

        #Drilling Cost Model
        self.Drilling_TimeMobilize_h = 336.0 #hr
        self.Drilling_CostMobilize_USD = 5130000.0 #USD
        self.Drilling_TimeRigWalk_h = 48.0 #hr
        self.Drilling_CostRigWalk_USD = 250000.0 #USD
        self.Drilling_TimeCut_hpm4 = 2.3e-11 #h/m4
        self.Drilling_TimeCut_hpm3 = 0.5 #h/m3
        self.Drilling_CostCut_USDph = 5082.1 #$/hr
        self.Drilling_CostCut_USDpm4 = 4.5e-5 #$/m4
        self.Drilling_CostCut_USDpm3 = 1000.0 #$/m3
        self.Drilling_TimeCasing_hpm4 = 2.0e-12 #h/m4
        self.Drilling_TimeCasing_hpm3 = 0.15 #h/m3
        self.Drilling_CostCasing_USDph = 4001.3 #$/hr
        self.Drilling_CostCasing_USDpm4 = 1.4 #$/m4
        self.Drilling_CostCasing_USDpm3 = 1460.0 #$/m3
        self.Drilling_Learning_ratio = 0.20 #learning rate
        
        #Stimulation Cost Model
        self.Fracking_Fixed_USD = 1233149.0 #$
        self.Fracking_Depth_USD = 638.0 #$/m
        self.Fracking_Hourly_USDph = 5973.3 #$/hr
        self.Fracking_Sand_USDpkg = 0.662 #$/kg #ISP
        self.Fracking_Fuel_USDpkWh = 0.451 #$/kWh
        self.Fracking_Water_USDpm3 = 22.2 #$/m3
        self.Fracking_PressureFactor_ratio = 0.8
        self.Fracking_PressureFactor_scale = 1e-8 #Pa
        self.Fracking_EquipmentWear_USDpm3 = 9.85 #$/m3
        self.Fracking_TimeFactor_hph = 1.95
        self.Fracking_StringRadius_m = 0.075 #m
        
        #Economic Parameters
        self.Cost_Planning_USD = -3.5e6 #USD, survey & planning costs
        self.Cost_Gridconnect_USD = -4.0e6 #USD, grid connection costs
        self.Cost_EquipmentCost_USDpkW = 2500.0 #%/kWe
        self.Cost_Grants_USD = 0.0e6 #USD, government grants
        self.Cost_Operations_ratio = 0.02 #ratio, opex vs capex
        self.Cost_Seismic_USDpMw = 2e-4 #$/Mw for $300M Mw 5.5 quake Pohang (Westaway, 2021) & $17.2B Mw 6.3 quake Christchurch (Swiss Re)
        self.Cost_Seismic_exp = 5.0 #$/Mw for $300M Mw 5.5 quake Pohang (Westaway, 2021) & $17.2B Mw 6.3 quake Christchurch (Swiss Re)
        self.Sales_ElectricWholesale_USDpkWh = 54.645e-3 #USD, market rate sales for electricity
        self.Sales_ElectricContract_USDpkWh = 165.590e-3 #USD, contract sales rate for electricity
        self.Sales_ElectricPurchase_USDpkWh = 291.300e-3 #USD, typical purchase price for customers
        self.Sales_HeatWholesale_USDpkWh = 10.000e-3 #USD, market rate sales for heat
        self.Sales_HeatContract_USDpkWh = 10.000e-3 #USD, contract rate sales for heat
        self.Sales_CapacityFactor_ratio = 0.95 #ratio, capacity factor for plant uptime
        self.Sales_ContractTime_yr = 15 #yr, length of contract for heat and electricity
        self.Equity_CostOfCapital_ratio = 0.12 #ratio, time value of money versus alternative investments
        self.Equity_Contribution_ratio = 0.40 #ratio, equity contribution assuming debt service reserve account funded by debt 
        self.Equity_CompanyShare_ratio = 0.02 #ratio, share of equity held by primary company
        self.Equity_SeriesAShare_ratio = 0.08 #ratio, share of equity held by series A investors
        self.Equity_SeriesBShare_ratio = 0.90 #ratio, share of equity held by series B investors
        self.Equity_BrokerFee_ratio = 0.050 #ratio, fee to setup project funds
        self.Debt_InterestRate_ratio = 0.080 #ratio, interest rate on debt and loans 
        self.Debt_ArrangementFee_ratio = 0.015 #ratio, fee to setup debt accounts
        self.Debt_ServiceAccount_ratio = 0.50 #principal and interest fraction of a year
        self.Debt_Tenor_yr = 30.0 #years, loan duration until maturity
        self.Tax_Inflation_ratio = 0.02 #inflation rate
        self.Tax_Corporate_ratio = 0.40 #ratio, tax rate on profits after losses 
        self.Tax_ProductionCredit_USDpkWh = 0.0 #0.015 #USD, tax credit for energy production
        self.Tax_CapitalCredit_ratio = 0.0  #0.40 #USD, tax credit for capital construction
        self.Tax_DirectPay_ratio = 0.0 #1.0 #ratio, ratio of tax credits that are elligble for direct pay (e.g., Inflation Reduction Act)
        self.Tax_DepreciationRate_ratio = 0.20
        self.CO2_Electric_kgpkWh = 0.22499 #kgCO2/kWhe, CO2 emitted from 'grid' power
        self.CO2_Heat_kgpkWh = 0.18293 #kgCO2/kWht, CO2 emitted from natural gas for heat 
        self.CO2_Geothermal_kgpkWh = 0.00154 #kgCO2/kWht, CO2 emitted from natural gas for heat
        
        #*******************************************************************************************************************************************************
        #relics for removal #!!! remove these before release
        self.Economics_ElectricitySales_USDpkWh = 0.1372 #$/kWh - customer electricity retail price 
        self.Economics_DrillingCost_USDpm = 2763.06 #$/m - Lowry et al, 2017 large diameter well baseline 
        self.Economics_AccessCost_USD = 590e3 #$ Lowry et al, 2017 large diameter well baseline 
        self.Economics_EquipmentCost_USDpkW = 2025.65 #$/kWe simplified from GETEM model 
        self.Economics_ExplorationCost_USDpm = 2683.41 #$/m simplified from GETEM model 
        self.Economics_Maintenance_USDpkWh = 0.03648 #$/kWh simplified from GETEM model 
        self.Economics_Seismic_USDpMw = 2e-4 #$/Mw for $300M Mw 5.5 quake Pohang (Westaway, 2021) & $17.2B Mw 6.3 quake Christchurch (Swiss Re) 
        self.Economics_Seismic_exp = 5.0 #$/Mw for $300M Mw 5.5 quake Pohang (Westaway, 2021) & $17.2B Mw 6.3 quake Christchurch (Swiss Re) 
        self.Well_HydraulicRadius_m = 0.0254*5.0 #m 
        self.Well_CasingRadius_m = 0.0254*5.5 # m 
        self.Well_BoreRadius_m = 0.0254*6.5 # m 
        self.StressMin_Coefficient_ratio = 0.5 #Pa/Pa 
        self.StressMax_Coefficient_ratio = 0.75 #Pa/Pa 
        self.StressMin_Azimuth_rad = 90.0*deg #rad 
        self.StressMin_Dip_rad = 0.0*deg #rad 
        self.StressMin_Variance_rad = 0.5*deg #rad
        self.Fracture_CompressibilityMin_1pPa = 2.0e-9 #1/Pa
        self.Fracture_CompressibilityNom_1pPa = 2.9e-8 #1/Pa
        self.Fracture_CompressibilityMax_1pPa = 10.0e-8 #1/Pa
        self.Proppant_CompressibilityMin_1pPa = 2.0e-9 #1/Pa
        self.Proppant_CompressibilityNom_1pPa = 2.9e-8 #1/Pa
        self.Proppant_CompressibilityMax_1pPa = 10.0e-8 #1/Pa
        self.Proppant_InitialPermeabilityMin_m2 = 10.0*darcy #m2
        self.Proppant_InitialPermeabilityNom_m2 = 100.0*darcy #m2
        self.Proppant_InitialPermeabilityMax_m2 = 300.0*darcy #m2
        self.Fracture_TortuosityMin_ratio = 0.80 #m/m
        self.Fracture_TortuosityNom_ratio = 0.90 #m/m
        self.Fracture_TortuosityMax_ratio = 1.00 #m/m
        self.Fluid_NominalDensity_kgpm3 = 965.0 #kg/m3
        self.Fluid_Viscosity_Pas = 0.2*cP #Pa-s
        self.Fracture_SlipLengthMin_ratio = 10.0**-3.0 #m/m
        self.Fracture_SlipLengthNom_ratio = 10.0**-2.0 #m/m
        self.Fracture_SlipLengthMax_ratio = 10.0**-1.2 #m/m
        self.Fracture_DilationSlipMin_ratio = 0.000 #m/m
        self.Fracture_DilationSlipNom_ratio = 0.200 #m/m
        self.Fracture_DilationSlipMax_ratio = 0.800 #m/m
        self.Fracture_HydraulicDilationMin_ratio = 0.0 #m/m
        self.Fracture_HydraulicDilationNom_ratio = 0.6 #m/m
        self.Fracture_HydraulicDilationMax_ratio = 2.0 #m/m
        self.Fracture_InitialHydraulicApertureMin_m = 0.00005 #m
        self.Fracture_InitialHydraulicApertureNom_m = 0.00010 #m 
        self.Fracture_InitialHydraulicApertureMax_m = 0.00020 #m 
        self.Shear_FrictionAngleMin_rad = 20.0*deg #rad
        self.Shear_FrictionAngleNom_rad = 35.0*deg #rad
        self.Shear_FrictionAngleMax_rad = 45.0*deg #rad
        self.Shear_CohesionMin_Pa = 1.0*1e6 #Pa
        self.Shear_CohesionNom_Pa = 7.0*1e6 #Pa
        self.Shear_CohesionMax_Pa = 15.0*1e6 #Pa
        self.Domain_StimulationLimit_steps = 1
        self.Strategy_SourceDirectory_path = os.getcwd()
        self.Strategy_TargetDirectory_path = os.getcwd()
        self.Strategy_Completion_type = 'DC' #type of completion; 'DC' = direct contact 'FE' = fracture exchange
        #*******************************************************************************************************************************************************
        
        self.re_init()
    
    def re_init(self):
        #shear modulus
        self.Rock_ShearModulus_Pa = self.Rock_YoungsModulus_Pa/(2.0*(1.0+self.Rock_PoissonRatio_ratio))
        
        #depth scaling
        if self.Strategy_Target_type.upper() == 'TEMP':
            self.Well_TargetDepth_m = (self.Well_TargetTemp_K-self.Air_Temperature_K)/self.Rock_ThermalGradient_Kpm
        else:
            self.Well_TargetTemp_K = self.Well_TargetDepth_m*self.Rock_ThermalGradient_Kpm+self.Air_Temperature_K #K
        self.Fluid_InjectionDensity_kgpm3 = 1.0/water.v_from_PT(P=self.Circulation_Backpressure_Pa,T=self.Circulation_InjectionTemperature_K)
        self.Fluid_InjectionViscosity_Pas = water.mu_from_T(T=self.Circulation_InjectionTemperature_K)
        self.Fluid_InjectionStatic_Pa = self.Fluid_InjectionDensity_kgpm3*g*self.Well_TargetDepth_m
        self.Fluid_ReservoirDensity_kgpm3 = 1.0/water.v_from_PT(P=self.Fluid_InjectionStatic_Pa,T=self.Circulation_InjectionTemperature_K)
        self.Fluid_ReservoirViscosity_Pas = water.mu_from_T(T=self.Well_TargetTemp_K)
        self.Fluid_ReservoirStatic_Pa = 0.5*(self.Fluid_InjectionDensity_kgpm3+self.Fluid_ReservoirDensity_kgpm3)*g*self.Well_TargetDepth_m
        if (self.Fluid_InjectionDensity_kgpm3<500) or (self.Fluid_ReservoirDensity_kgpm3<500):
            print("NOTICE: injected water is not a liquid at surface or at bottom-hole")
        
        #intervals, clusters, and well count (per unit of pattern)
        if self.Strategy_Design_type.upper() == 'AGS':
            self.Well_Intervals_count = 1
            self.Well_Clusters_count = 1
            self.Well_InjectorCount_wells = self.Well_ProducerCount_wells
        elif self.Strategy_Design_type.upper() == 'FGS':
            self.Well_Intervals_count = max([1,int(self.Stimulation_Clusters_clusters)])
            self.Well_Clusters_count = 1
            self.Well_InjectorCount_wells = 1
        elif self.Strategy_Design_type.upper() == 'CGS':
            self.Well_Intervals_count = max([1,int(self.Stimulation_Clusters_clusters)])
            self.Well_Clusters_count = 1
            self.Well_InjectorCount_wells = int(self.Well_Pattern_count*max([1,int(self.Well_ProducerCount_wells-1)]))
        elif self.Strategy_Design_type.upper() == 'O&G':
            self.Well_Intervals_count = 1
            self.Well_Clusters_count = max([1,int(self.Stimulation_Clusters_clusters)])
            self.Well_InjectorCount_wells = max([1,int(self.Well_ProducerCount_wells-1)])
        else: # 'EGS' 'DCM'
            self.Well_Intervals_count = 1
            self.Well_Clusters_count = max([1,int(self.Stimulation_Clusters_clusters)])
            self.Well_InjectorCount_wells = 1
        
        #well scaling
        self.Well_ProductionRadius_m = 0.5*self.Well_ProductionDiameter_m
        self.Well_SegmentCount_units = max([1,int(self.Well_SegmentCount_units)])
        
        #stress conditions
        s123 = np.zeros(3,dtype=float)
        s123[0] = self.Rock_Density_kgpm3*g*self.Well_TargetDepth_m
        s123[2] = self.Stress_ShminToSv_ratio*(s123[0]-self.Fluid_ReservoirStatic_Pa*self.Rock_Biot_ratio) + self.Fluid_ReservoirStatic_Pa*self.Rock_Biot_ratio
        s123[1] = self.Stress_ShmaxToShmin_ratio*(s123[2]-self.Fluid_ReservoirStatic_Pa*self.Rock_Biot_ratio) + self.Fluid_ReservoirStatic_Pa*self.Rock_Biot_ratio
        tensor = stress.cauchy()
        tensor.set_sigG_from_Principal(s123[2], s123[1], s123[0], self.Stress_ShminAzimuth_rad, self.Stress_ShminDip_rad)
        str_dip = tensor.get_conjugates(plots=True)
        s123.sort(axis=0)
        self.Stress_S1_Pa = s123[2]
        self.Stress_S2_Pa = s123[1]
        self.Stress_S3_Pa = s123[0]
        self.JointSets_Set1Strike_rad = str_dip[0,0]
        self.JointSets_Set1Dip_rad = str_dip[0,1]
        self.JointSets_Set2Strike_rad = str_dip[1,0]
        self.JointSets_Set2Dip_rad = str_dip[1,1]
        self.JointSets_Set3Strike_rad = str_dip[2,0]
        self.JointSets_Set3Dip_rad = str_dip[2,1]
        
        #simplistic calculations
        self.Well_DrilledLength_m = self.Well_TargetDepth_m + self.Well_ProducerProportion_ratio*self.Well_DeviatedLength_m*(1 - 0.5*np.sin(self.Well_Dip_rad)) 
        self.Stimulation_Pressure_Pa = ((10.7/0.9e-3)*(self.Well_DrilledLength_m*self.Fluid_InjectionViscosity_Pas*self.Fluid_InjectionDensity_kgpm3*g*self.Stimulation_InjectionRate_m3ps**1.852)
                                        /(self.Well_Roughness_metric**1.852*(self.Well_ProductionDiameter_m-2.0*self.Well_WallThickness_m)**4.87) 
                                        +self.Stress_S3_Pa + self.Tension_Cohesion_Pa - self.Fluid_InjectionStatic_Pa)
        
#standard outputs
class payoff:
    def __init__(self):
        #keymost parameters
        self.Profit_ProjectIRR_ratio = 0.0
        self.Profit_LCOE_USDpMWh = 0.0
        self.Profit_LCOH_USDpMWh = 0.0
        self.Profit_NetPresentValue_USD = 0.0
        self.Power_RunningNet_kW = 0.0
        self.Power_RunningGross_kW = 0.0
        self.Power_RunningThermal_kW = 0.0
        self.Power_RunningPump_kW = 0.0
        self.Power_ProductiveLife_yr = 0.0
        
        #temperatures
        self.Temperature_Injection_C = 0.0
        self.Temperature_Reservoir_C = 0.0
        self.Temperature_Decline_Cpyr = 0.0
        
        #flow
        self.Pressure_Stimulation_MPa = 0.0
        self.Pressure_Injection_MPa = 0.0
        self.Pressure_Production_MPa = 0.0
        self.Rate_Injection_kgps = 0.0
        self.Rate_Production_kgps = 0.0
        self.Rate_Leakoff_kgps = 0.0
        
        #economics
        self.Sales_Heat_USD = 0.0
        self.Sales_Electric_USD = 0.0
        self.Cost_Pumping_USD = 0.0
        self.Cost_Facilities_USD = 0.0 #USD, power plant cost
        self.Cost_Fracking_USD = 0.0 #USD, cost to stimulate
        self.Cost_Drilling_USD = 0.0 #USD, cost of drilling
        self.Cost_Time_yr = 0.0 #yr, time until first power production
        self.Profit_EquityIRR_ratio = 0.0
        self.Profit_SeriesAIRR_ratio = 0.0
        self.Profit_SeriesBIRR_ratio = 0.0
        
        #seimsicity
        self.Stress_s3_MPa = 0.0
        self.Stress_s2_MPa = 0.0
        self.Stress_s1_MPa = 0.0
        self.Stress_pp_MPa = 0.0
        self.Seismicity_MaxQuake_Mw = 0.0
        self.Seismicity_NumQuake_ea = 0.0
        self.Fractures_ShearStim_ea = 0.0
        self.Fractures_HydroStim_ea = 0.0
        self.Fractures_Hydroprop_ea = 0.0        
        
        #tidbits
        self.Drilling_TotalLength_m = 0.0
        self.Drilling_CementVolume_m3 = 0.0
        self.Drilling_CasingWeight_kgpm = 0.0
        self.Stimulation_SandMass_kg = 0.0
        self.Stimulation_WaterVolume_m3 = 0.0
        self.Circulation_WorkingVolume_m3 = 0.0
        self.Circulation_Recovery_pct = 0.0
        self.Circulation_MakeupWater_m3 = 0.0
        self.CO2_Offset_MtCO2eq = 0.0
        
        
        
        
        
    #     #Emissions offsetting
    #     1e-9*np.sum(self.Power_Thermal_kWh*(self.CO2_Heat_kgpkWh-self.CO2_Geothermal_kgpkWh)*keep 
    #                                 +(self.Power_Gross_kWh+self.Power_Pumping_kWh)*(self.CO2_Electric_kgpkWh-self.CO2_Geothermal_kgpkWh)*keep)
        
                
    #         y = data['Injection_Rate_m3ps']; yname = 'Circulation Rate (m3/s)'
    # # y = data['Rock_ThermalConductivity_WpmK']; yname = 'Thermal Conductivity (W/m-K)'
    # # y = data['Stimulation_Clusters_clusters']; yname = 'Number of Fractures (stages)' ; ylog = False
    
    # # x = data['Well_Spacing_m']; name = 'Well Spacing (m)' ; xlog = False
    # # x = data['Stimulation_ProppantConcentration_m3pm3']; name = 'Proppant Concentration (m3/m3)' ; xlog = False
    # # x = data['Proppant_InitialPermeabilityNom_m2']; name = 'Proppant Permeability (m2)' ; xlog = True
    # # x = data['Stimulation_Clusters_clusters']; name = 'Number of Fractures (stages)' ; xlog = False
    # # x = data['Well_DeviatedLength_m']; name = 'Lateral Length (m)' ; xlog = False
    # # x = data['Rock_ThermalGradient_Kpm']*1e3; name = 'Thermal Gradient (C/km)' ; xlog = False
    # # x = data['Stimulation_InjectionRate_m3ps']; name = 'Stimulation Rate (m3/s)' ; xlog = True
    # # x = data['Well_BoreRadius_m']; name = 'Borehole Radius (m)' ; xlog = False
    # # x = data['StressMin_Coefficient_ratio']; name = 'Stress Condition (0.5-0.8 NF; 0.8-1.2 SS or RF)' ; xlog = False
    # # producers = data['Well_ProducerCount_wells']
    # # injectors = np.max([producers-1,np.ones(len(data))],axis=0)
    # # x = injectors+producers; name = 'Number of Wells (Injectors & Producers)' ; xlog = False
    # # x = data['Well_TargetDepth_m']; name = 'Depth (m)' ; xlog = False #l3
    # # x = data['Stimulation_TargetRadius_m']; name = 'Fracture Radius (m)' ; xlog = False #l3
    # x = data['Injection_BottomPressure_MPa']-data['Production_BottomPressure_MPa']; name = 'Fault Pressure Loss (MPa)'; xlog = True
    # x = data['Fracture_Conductivity_m2m']*3.324e15; name = 'Fracture Conductivity (mD-ft)'; xlog = True
        
        
        
        