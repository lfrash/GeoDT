# -*- coding: utf-8 -*-
""" 
Classes for calculating cashflow for a geothermal development project:
    1. Debt, equity, fees, and taxes
    2. IRR, LCOE, and NPV
    3. Scheduled power and heat delivery
    4. Two phase pricing schedule (introductory & long term)
    5. Attributes include annual cash flows
    6. Well cost
    7. Drilling cost
"""
#imports
import numpy as np
if __package__ is None or __package__ == '':
    import iogt
    import wells
    import properties
    from units import *
else:
    from . import iogt
    from . import wells
    from . import properties
    from .units import *
import pylab
water = properties.water()

#units
deg=deg
g=g
yr=yr

#**********************************************************************************************************
### economics functions
#**********************************************************************************************************
class cashflow():
    def __init__(self): #defaults to "Durham EGS 3P 3I lpf v2.xlsx - CF GBP"
        #project parameters
        self.Cost_PlanningTime_yr = 0.99 #yr, final investment decision
        self.Cost_Capital_yr = 2.99 #yr, construction time including planning
        self.Cost_Lifespan_yr = 32 #yr, total project life
        self.Cost_Planning_USD = -3.5430875e6 #USD, survey & planning costs
        self.Cost_Capital_USD = -225.025460e6 #USD, drilling, stimulation, grid, facilities, and grants
        self.Cost_Operations_ratio = 0.02 #ratio, opex vs capex
        #sales pricing
        self.Sales_ElectricWholesale_USDpkWh = 54.645e-3 #USD, market rate sales for electricity
        self.Sales_ElectricContract_USDpkWh = 165.590e-3 #USD, contract sales rate for electricity
        self.Sales_ElectricPurchase_USDpkWh = 291.300e-3 #USD, typical purchase price for customers
        self.Sales_HeatWholesale_USDpkWh = 10.000e-3 #USD, market rate sales for heat
        self.Sales_HeatContract_USDpkWh = 10.000e-3 #USD, contract rate sales for heat
        self.Sales_CapacityFactor_ratio = 0.95 #ratio, capacity factor for plant uptime
        self.Sales_ContractTime_yr = 15 #yr, length of contract for heat and electricity
        #equity, debt, and taxes
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
        #performance, noting that length of these series is shorter than Capex array
        self.Power_Gross_kWh = np.asarray([0,0,0,311220.7353,292546.4382,277972.1724,260977.3087,247599.0859,232420.428,219786.6588,
                              208207.0261,196763.7449,185459.8462,174297.3886,162840.5088,151213.9588,139890.924,
                              128853.0469,118084.7095,107572.5474,97305.06136,86414.00885,75732.30733,65867.86571,
                              56509.57654,47513.06953,38802.21546,30333.07111,22078.4974,14020.82846,6148.07003,
                              -1548.207607,-9074.102275,-16434.13404,-23631.72431])*1e3
        self.Power_Thermal_kWh = np.asarray([0,0,0,0,24382.75325,32510.33766,56893.09091,65020.67532,91689.31169,100578.8571,
                                  100578.8571,100578.8571,100578.8571,100578.8571,100578.8571,100578.8571,100578.8571,
                                  100578.8571,100578.8571,100578.8571,100578.8571,100578.8571,100578.8571,100578.8571,
                                  100578.8571,100578.8571,100578.8571,100578.8571,100578.8571,100578.8571,100578.8571,
                                  100578.8571,100578.8571,100578.8571,100578.8571])*1e3
        self.Power_Pumping_kWh = np.asarray([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])*1e3
        self.Power_Heat_kWh = self.Power_Gross_kWh/0.15 + self.Power_Thermal_kWh
        #CO2 offset, https://www.gov.uk/government/publications/greenhouse-gas-reporting-conversion-factors-2024
        self.CO2_Electric_kgpkWh = 0.22499 #kgCO2/kWhe, CO2 emitted from 'grid' power
        self.CO2_Heat_kgpkWh = 0.18293 #kgCO2/kWht, CO2 emitted from natural gas for heat 
        self.CO2_Geothermal_kgpkWh = 0.00154
        #exchange rate
        self.Currency_Exchange_USDpGBP = 1.3 #USD, currency exchange rate
        #calculated parameters
        self.re_init()
    def re_init(self):
        #refernce tables
        span = int(self.Cost_Lifespan_yr+1+int(self.Cost_Capital_yr))
        self.Equity_WACC_ratio = self.Equity_Contribution_ratio*self.Equity_CostOfCapital_ratio+(
            (1.0-self.Equity_Contribution_ratio)*self.Debt_InterestRate_ratio*(1.0-self.Tax_Corporate_ratio))
        self.Inflation_Table_series = (1.0+self.Tax_Inflation_ratio)**np.linspace(0,span-1,span)
        self.Discount_Table_series = 1.0/((1.0+self.Equity_WACC_ratio)**np.linspace(0,span-1,span))
        #series
        self.Cost_Capex_USD = np.zeros(span)
        self.Cost_Opex_USD = np.zeros(span)
        a = 1+int(self.Cost_PlanningTime_yr-0.001)
        b = 1+int(self.Cost_Capital_yr-0.001)
        self.Cost_Capex_USD[a:b] += ((self.Cost_Capital_USD)/(b-a))*np.ones(b-a)
        for i in range(b,span):
            self.Cost_Opex_USD[i] = np.sum(self.Cost_Capex_USD[:i])*self.Cost_Operations_ratio
        self.Cost_Capex_USD[0:a] += (self.Cost_Planning_USD/a)*np.asarray(a)
        #capital and operations costs with inflation
        self.Cost_Capex_USD = self.Cost_Capex_USD*self.Inflation_Table_series
        self.Cost_Opex_USD = self.Cost_Opex_USD*self.Inflation_Table_series
        #depreciation table
        self.Tax_Depreciation_USD = np.zeros(span)
        for i in range(0,span):
            Depreciation_StraightLine_yr = np.max([1,int(1/self.Tax_DepreciationRate_ratio)]) 
            capital = np.ones(Depreciation_StraightLine_yr)*self.Cost_Capex_USD[i]/Depreciation_StraightLine_yr
            j = np.min([i+Depreciation_StraightLine_yr,span])
            self.Tax_Depreciation_USD[i:j] += capital[:j-i]
        #revenue from heat and electricity sales
        j = np.min([span,b+int(self.Sales_ContractTime_yr)])
        self.Sales_Electric_USDpkWh = self.Sales_ElectricWholesale_USDpkWh*self.Inflation_Table_series
        self.Sales_Electric_USDpkWh[:j] = self.Sales_ElectricContract_USDpkWh*self.Inflation_Table_series[:j]
        self.Sales_Heat_USDpkWh = self.Sales_HeatWholesale_USDpkWh*self.Inflation_Table_series
        self.Sales_Heat_USDpkWh[:j] = self.Sales_HeatContract_USDpkWh*self.Inflation_Table_series[:j]
        self.Sales_Energy_USD = self.Power_Gross_kWh*self.Sales_Electric_USDpkWh + self.Power_Thermal_kWh*self.Sales_Heat_USDpkWh
        #pumping losses with electric power and production tax credits
        buy = np.argmin(np.asarray([self.Sales_Electric_USDpkWh,self.Sales_ElectricPurchase_USDpkWh*np.ones(span)]),axis=0)
        self.Cost_Electric_USDpkWh = (self.Sales_ElectricPurchase_USDpkWh*buy+self.Sales_Electric_USDpkWh*(1.0-buy))*self.Inflation_Table_series
        self.Tax_ProductionCredit_USD = (self.Power_Gross_kWh+self.Power_Pumping_kWh*(1.0-buy))*self.Tax_ProductionCredit_USDpkWh
        self.Cost_Energy_USD = self.Power_Pumping_kWh*self.Cost_Electric_USDpkWh
        self.Tax_ProductionCredit_USD[self.Tax_ProductionCredit_USD<0] = 0.0
        self.Tax_CapitalCredit_USD = -1.0*self.Cost_Capex_USD*self.Tax_CapitalCredit_ratio
        #cumulative tax credits (direct pay and liability)
        self.Tax_CumulativeCredit_USD = np.zeros(span)
        self.Tax_DirectPay_USD = np.zeros(span)
        self.Tax_CumulativeCredit_USD[0] = self.Tax_ProductionCredit_USD[0]+(self.Tax_CapitalCredit_USD[0])*(1.0-self.Tax_DirectPay_ratio)
        self.Tax_DirectPay_USD[0] = (self.Tax_CapitalCredit_USD[0])*self.Tax_DirectPay_ratio*(1.0-self.Equity_WACC_ratio)**b
        for i in range(0,span):
            self.Tax_CumulativeCredit_USD[i] = self.Tax_ProductionCredit_USD[i]+(self.Tax_CapitalCredit_USD[i])*(1.0-self.Tax_DirectPay_ratio)+self.Tax_CumulativeCredit_USD[i-1]
            self.Tax_DirectPay_USD[i] = (self.Tax_CapitalCredit_USD[i])*self.Tax_DirectPay_ratio*(1.0-self.Equity_WACC_ratio)**np.min([1,b-i])
        #profit before tax with equity, debt, servicing, interest, and fees
        self.Debt_Equity_USD = np.zeros(span)
        self.Debt_CashFlow_USD = np.zeros(span)
        self.Debt_Cumulative_USD = np.zeros(span)
        self.Debt_Interest_USD = np.zeros(span)
        self.Debt_Reserve_USD = np.zeros(span)
        self.Debt_Arrangement_USD = np.zeros(span)
        self.Profit_CashFlow_USD = np.zeros(span)
        self.Tax_LossCarryForward_USD = np.zeros(span)
        K = self.Cost_Capex_USD+self.Cost_Opex_USD+self.Sales_Energy_USD+self.Cost_Energy_USD+self.Tax_DirectPay_USD
        rB = self.Equity_BrokerFee_ratio
        rA = self.Debt_ArrangementFee_ratio
        rE = self.Equity_Contribution_ratio
        rR = self.Debt_ServiceAccount_ratio
        rI = self.Debt_InterestRate_ratio
        T = self.Debt_Tenor_yr
        x1 = rE*(1.0-rB)
        x2 = 1.0-rB*rE
        self.Debt_CashFlow_USD[0] = 0.0
        self.Debt_Cumulative_USD[0] = 0.0
        self.Debt_Equity_USD[0] = K[0]/(1-rB)
        self.Debt_Interest_USD[0] = 0.0
        self.Debt_Reserve_USD[0] = 0.0
        self.Debt_Arrangement_USD[0] = self.Debt_Equity_USD[0]*rB
        self.Profit_CashFlow_USD[0] = K[0]+self.Debt_Interest_USD[0]+self.Debt_Arrangement_USD[0]
        self.Tax_LossCarryForward_USD[0] = self.Profit_CashFlow_USD[0]-self.Cost_Capex_USD[0]+self.Tax_Depreciation_USD[0]
        for i in range(1,span):
            #debt calculations
            Ci = self.Debt_Cumulative_USD[i-1]
            Ri = self.Debt_Reserve_USD[i-1]
            #development phase
            if i < b:
                self.Debt_CashFlow_USD[i] = ((K[i]+Ci*rI+Ci*rI*rR-Ri)*x2 - (K[i]+Ci*rI)*x1)/(rI*x1+rA*x1-rI*x2-rA*x2+x2-rI*rR*x2)
                self.Debt_Cumulative_USD[i] = self.Debt_CashFlow_USD[i]+Ci
                self.Debt_Equity_USD[i] = (K[i]+self.Debt_Cumulative_USD[i]*rI+self.Debt_CashFlow_USD[i]*rA)*rE/x2
                self.Debt_Interest_USD[i] = self.Debt_Cumulative_USD[i]*rI 
                self.Debt_Reserve_USD[i] = self.Debt_Interest_USD[i]*rR
            #production phase assuming positive cash flow
            else:
                #debt calculations
                self.Debt_Equity_USD[i] = 0.0
                self.Debt_CashFlow_USD[i] = -1.0*np.min(self.Debt_Cumulative_USD[:i])/T
                self.Debt_Cumulative_USD[i] = self.Debt_CashFlow_USD[i]+Ci
                self.Debt_Interest_USD[i] = self.Debt_Cumulative_USD[i]*rI 
                self.Debt_Reserve_USD[i] = (self.Debt_Interest_USD[i]-self.Debt_CashFlow_USD[i])*rR
                self.Debt_Arrangement_USD[i] = np.min([0.0,self.Debt_Equity_USD[i]])*rB+np.min([0.0,self.Debt_CashFlow_USD[i]])*rA
                self.Profit_CashFlow_USD[i] = K[i]+self.Debt_Interest_USD[i]+self.Debt_Arrangement_USD[i]
                #payoff loan early if beginning to take on losses other than through capital investment
                if (self.Profit_CashFlow_USD[i-1]<0.0) and (self.Debt_Cumulative_USD[i-1]<-1.0) and (self.Cost_Capex_USD[i-1]>-1.0):
                    self.Debt_Equity_USD[i] = self.Debt_Cumulative_USD[i-1] #-self.Debt_Reserve_USD[i-1]
                    self.Debt_CashFlow_USD[i] = -1.0*self.Debt_Equity_USD[i]
                    self.Debt_Cumulative_USD[i] = self.Debt_CashFlow_USD[i]+Ci
                #debt is paid at end of tenor
                elif self.Debt_Cumulative_USD[i-1]>-1.0:
                    self.Debt_Equity_USD[i] = 0.0
                    self.Debt_CashFlow_USD[i] = 0.0
                    self.Debt_Cumulative_USD[i] = 0.0
                #production phase if negative cash flow (assumes debt not yet paid)
                elif (self.Profit_CashFlow_USD[i] < 0):
                    self.Debt_Equity_USD[i] = 0.0
                    self.Debt_CashFlow_USD[i] = (K[i]-np.min(self.Debt_Cumulative_USD[:i])/T+Ci*rI+Ci*rI*rR)/(1.0-rI-rA-rI*rR)
                    self.Debt_Cumulative_USD[i] = self.Debt_CashFlow_USD[i]+Ci
                self.Debt_Interest_USD[i] = self.Debt_Cumulative_USD[i]*rI
                if self.Debt_Cumulative_USD[i] < -1.0:
                    self.Debt_Reserve_USD[i] = (self.Debt_Interest_USD[i]-self.Debt_CashFlow_USD[i])*rR
                else:
                    self.Debt_Reserve_USD[i] = self.Debt_Interest_USD[i]*rR
            #cash flow
            self.Debt_Arrangement_USD[i] = np.min([0.0,self.Debt_Equity_USD[i]])*rB+np.min([0.0,self.Debt_CashFlow_USD[i]])*rA
            self.Profit_CashFlow_USD[i] = K[i]+self.Debt_Interest_USD[i]+self.Debt_Arrangement_USD[i]
            self.Tax_LossCarryForward_USD[i] = self.Profit_CashFlow_USD[i]-self.Cost_Capex_USD[i]+self.Tax_Depreciation_USD[i]+np.min([0.0,self.Tax_LossCarryForward_USD[i-1]])
            
        #taxes
        self.Tax_Levied_USD = np.zeros(span)
        self.Profit_BeforeTax_USD = self.Profit_CashFlow_USD-self.Cost_Capex_USD
        self.Tax_LossCarryForward_USD[self.Tax_LossCarryForward_USD>0] = 0.0
        self.Tax_Levied_USD[0] = -(self.Profit_BeforeTax_USD[0]+self.Tax_Depreciation_USD[0])*self.Tax_Corporate_ratio
        self.Tax_Levied_USD[1:] = -(self.Profit_BeforeTax_USD[1:]+self.Tax_Depreciation_USD[1:]+self.Tax_LossCarryForward_USD[:-1])*self.Tax_Corporate_ratio
        self.Tax_Levied_USD[self.Tax_Levied_USD>0] = 0.0
        self.Profit_AfterTax_USD = np.zeros(span)
        for i in range(0,span):
            tax = np.min([0.0,self.Tax_Levied_USD[i]+self.Tax_CumulativeCredit_USD[i]])
            self.Tax_CumulativeCredit_USD[i:] = self.Tax_CumulativeCredit_USD[i:]+tax
            self.Tax_CumulativeCredit_USD[self.Tax_CumulativeCredit_USD<0] = 0.0
            self.Tax_Levied_USD[i] = tax
            self.Profit_AfterTax_USD[i] = self.Profit_BeforeTax_USD[i]+self.Tax_Levied_USD[i]
        
        #cashflows
        self.Profit_Equity_USD = np.zeros(span)
        self.Profit_CashFlowTaxed_USD = self.Profit_CashFlow_USD+self.Tax_Levied_USD
        self.Profit_Equity_USD[0] = self.Profit_AfterTax_USD[0]+self.Cost_Capex_USD[0]-self.Debt_CashFlow_USD[0]+self.Debt_Reserve_USD[0]
        self.Profit_Equity_USD[1:] = self.Profit_AfterTax_USD[1:]+self.Cost_Capex_USD[1:]-self.Debt_CashFlow_USD[1:]+self.Debt_Reserve_USD[1:]-self.Debt_Reserve_USD[:-1]
        keep = np.zeros(span)
        keep[self.Debt_Equity_USD<-1.0] = 1
        keep[self.Profit_AfterTax_USD>0.0] = 1
        keep[self.Debt_Cumulative_USD<-0.1] = 1
        self.Profit_Equity_USD[keep<1] = 0.0
        self.Profit_Company_USD = self.Profit_Equity_USD*self.Equity_CompanyShare_ratio
        if self.Equity_SeriesAShare_ratio > 0:
            self.Profit_SeriesA_USD = self.Profit_Equity_USD*(1.0-self.Equity_CompanyShare_ratio)*self.Equity_SeriesAShare_ratio/(self.Equity_SeriesAShare_ratio+self.Equity_SeriesBShare_ratio)
        else:
            self.Profit_SeriesA_USD = np.zeros(span)
        if self.Equity_SeriesBShare_ratio > 0:
            self.Profit_SeriesB_USD = self.Profit_Equity_USD*(1.0-self.Equity_CompanyShare_ratio)*self.Equity_SeriesBShare_ratio/(self.Equity_SeriesAShare_ratio+self.Equity_SeriesBShare_ratio)
        else:
            self.Profit_SeriesB_USD = np.zeros(span)
        plug = np.ones(b)
        plug[0] = 0.0
        self.Profit_Company_USD[:b] = 0.0
        self.Profit_SeriesA_USD[:b] = (1-plug)*self.Debt_Equity_USD[:b]
        self.Profit_SeriesB_USD[:b] = (plug)*self.Debt_Equity_USD[:b]
        
        #summations
        self.keep = keep
        self.Sales_Heat_USD = np.sum(self.Power_Thermal_kWh*self.Sales_Heat_USDpkWh*keep)
        self.Sales_Electric_USD = np.sum(self.Power_Gross_kWh*self.Sales_Electric_USDpkWh*keep)
        self.Cost_Pumping_USD = np.sum(self.Cost_Energy_USD*keep)
        self.Cost_Operations_USD = np.sum(self.Cost_Opex_USD*keep)
        self.Cost_Credits_USD = np.sum(self.Tax_CumulativeCredit_USD*keep)
        self.Power_ProductiveLife_yr = np.sum(keep)
        self.Power_RunningGross_kW = np.sum(self.Power_Gross_kWh*keep/(365.25*24))/np.sum(keep)
        self.Power_RunningThermal_kW = np.sum(self.Power_Thermal_kWh*keep/(365.25*24))/np.sum(keep)
        self.Power_RunningPump_kW = np.sum(self.Power_Pumping_kWh*keep/(365.25*24))/np.sum(keep)
        self.Power_RunningNet_kW = np.sum((self.Power_Gross_kWh+self.Power_Pumping_kWh)*keep/(365.25*24))/np.sum(keep)
        
        #IRR and NPV
        self.Profit_EquityIRR_ratio = self.irr(self.Profit_Equity_USD)
        self.Profit_SeriesAIRR_ratio = self.irr(self.Profit_SeriesA_USD)
        self.Profit_SeriesBIRR_ratio = self.irr(self.Profit_SeriesB_USD)
        self.Profit_ProjectIRR_ratio = self.irr(self.Profit_CashFlowTaxed_USD[keep>0])
        self.Profit_NetPresentValue_USD = np.sum(self.Profit_CashFlowTaxed_USD[keep>0]/self.Discount_Table_series[keep>0])
        
        #LCOE and LCOH
        self.Profit_LCOE_USDpMWh = -1e3*(np.sum((self.Cost_Capex_USD+self.Cost_Opex_USD)*self.Discount_Table_series*keep)
                                         /np.max([1,np.sum((self.Power_Gross_kWh+self.Power_Pumping_kWh)*self.Discount_Table_series*keep)]))
        self.Profit_LCOH_USDpMWh = -1e3*(np.sum((self.Cost_Capex_USD+self.Cost_Opex_USD)*self.Discount_Table_series*keep)
                                        /np.max([1,np.sum(self.Power_Heat_kWh*self.Discount_Table_series*keep)]))
        
        #Emissions offsetting
        self.CO2_Offset_MtCO2eq = 1e-9*np.sum(self.Power_Thermal_kWh*(self.CO2_Heat_kgpkWh-self.CO2_Geothermal_kgpkWh)*keep 
                                    +(self.Power_Gross_kWh+self.Power_Pumping_kWh)*(self.CO2_Electric_kgpkWh-self.CO2_Geothermal_kgpkWh)*keep)
        
    def irr(self,cashflows):
        rate = 0.0
        lo = -10.0
        hi = 100.0
        for i in range(0,25):
            npv = (1 + rate)**np.linspace(0,len(cashflows)-1,len(cashflows))
            npv[npv==0] = 0.1
            try:
                npv = np.sum(cashflows/npv)
            except:
                print('error in irr calculation')
                print(npv)
                npv = 0.0
            if npv > 0.1:
                lo = rate
            else:
                hi = rate
            rate = (lo + hi)/2
            # print(rate)
        return rate

#**********************************************************************************************************
### calcuate drilling costs
#**********************************************************************************************************
# def drilling(s=iogt.setup(),w=wells.gen_wells(iogt.setup()),visuals=False):
def drilling(s=[],w=[],visuals=False):
    #well information extractor
    count = 0
    for i in range(0,len(w)):
        count = np.max([count,w[i].pID])
    wids = np.linspace(0,count,count+1) #array of well IDs
        
    #learning table
    learning = (wids+1)**(np.log(1.0-s.Drilling_Learning_ratio)/np.log(2.0))
    
    #running totals
    Drilling_Time_h = s.Drilling_TimeMobilize_h
    Drilling_Cost_USD = s.Drilling_CostMobilize_USD
    if visuals:
        cost_tracker = [[0.0,0.0]]
        cost_tracker += [[Drilling_Time_h,Drilling_Cost_USD]]
    
    #estimated cost for drilling campaign
    for i in range(0,len(w)):
        #segment midpoint depth
        middepth = s.Well_TargetDepth_m - (0.5*(w[i].c0 + w[i].c1))[2] #m
        cutvolume = w[i].leg*np.pi*w[i].rc**2 #m3
        
        #drilling time-cost
        # time_cut = ((s.Drilling_TimeCut_hpm4*middepth**3.0 + s.Drilling_TimeCut_hpm3) * cutvolume) * learning[w[i].pID]
        time_cut = ((s.Drilling_TimeCut_hpm2*middepth**2.0 + s.Drilling_TimeCut_hpm) * w[i].leg) * learning[w[i].pID]
        cost_cut = time_cut*s.Drilling_CostCut_USDph
        
        #drilling cut-cost
        # cost_cut += (s.Drilling_CostCut_USDpm4*middepth**2.0 + s.Drilling_CostCut_USDpm3) * cutvolume
        cost_cut += (s.Drilling_CostCut_USDpm4*middepth**1.0 + s.Drilling_CostCut_USDpm3) * cutvolume
        
        #total time
        Drilling_Time_h += time_cut
        Drilling_Cost_USD += cost_cut
        #visuals
        if visuals: 
            cost_tracker += [[Drilling_Time_h,Drilling_Cost_USD]]
        
        #check if cased
        time_case = 0.0
        cost_case = 0.0
        if not(w[i].type in ['screen']):
            #casing time
            # time_case = ((s.Drilling_TimeCasing_hpm4*middepth**3.0 + s.Drilling_TimeCasing_hpm3) * cutvolume) * learning[w[i].pID]
            time_case = s.Drilling_TimeCasing_hpm * w[i].leg * learning[w[i].pID]
            cost_case = time_cut*s.Drilling_CostCasing_USDph
            
            #casing cost
            cost_case += ((s.Drilling_CostCasing_USDpm4*middepth**1.0 + s.Drilling_CostCasing_USDpm3) * cutvolume)
            
        #total time
        Drilling_Time_h += time_case
        Drilling_Cost_USD += cost_case
        #visuals
        if visuals: 
            cost_tracker += [[Drilling_Time_h,Drilling_Cost_USD]]
    
    #final costs
    Drilling_Time_h += np.sum(wids*s.Drilling_TimeRigWalk_h)
    Drilling_Cost_USD += np.sum(wids*s.Drilling_CostRigWalk_USD)
    
    #visuals
    if visuals:
        cost_tracker += [[Drilling_Time_h,Drilling_Cost_USD]]
        cost_tracker = np.asarray(cost_tracker)
        fig = pylab.figure(figsize=(10.0,5.0),dpi=100)
        ax = fig.add_subplot(111)
        ax.plot(cost_tracker[:,0]/24,cost_tracker[:,1],'.-',color='blue')
        ax.set_ylabel('Drilling Cost (USD)',fontsize=10)
        ax.set_xlabel('Drilling Time (d)',fontsize=10)
    
    #result
    s.Cost_Drilling_h = Drilling_Time_h
    s.Cost_Drilling_USD = -Drilling_Cost_USD
    return Drilling_Time_h, -Drilling_Cost_USD

#**********************************************************************************************************
### calcuate stimulation costs
#**********************************************************************************************************
def fracking(s=[],w=[],visuals=False):
    #no injector
    if (s.Well_InjectorCount_wells == 0) or (s.Well_Clusters_count == 0):
        Stimulation_Time_h = 0.0
        Stimulation_Cost_USD = 0.0
    else:
        #common information
        Pstim = s.Stimulation_Pressure_Pa
        Qstim = s.Stimulation_InjectionRate_m3ps
        c_v = 0.0
        for i in range(0,len(w)):
            if w[i].type in ['perfcluster','procluster']:
                c_v += s.Well_Clusters_count * s.Stimulation_Volume_m3pfrac #slurry volume with leakoff
        c_m = c_v * s.Stimulation_ProppantConcentration_m3pm3 * s.Stimulation_ProppantDensity_kgpm3 #total proppant mass
        
        #fixed costs (mobilization, pad, parts)
        Stimulation_Cost_USD = s.Fracking_Fixed_USD
        
        #depth related fees
        Stimulation_Cost_USD += s.Fracking_Depth_USD*s.Well_DrilledLength_m 
        
        #volume costs
        Stimulation_Cost_USD += s.Fracking_Water_USDpm3*c_v
        Stimulation_Cost_USD += s.Fracking_Sand_USDpkg*c_m
        
        #time costs
        Stimulation_Time_h = s.Fracking_TimeFactor_hph*c_v/(Qstim*60*60)
        Stimulation_Cost_USD += s.Fracking_Fuel_USDpkWh*Stimulation_Time_h*Pstim*Qstim*1e-3
        Stimulation_Cost_USD += Stimulation_Time_h*s.Fracking_Hourly_USDph
        
        #pressure-volume costs 
        Stimulation_Cost_USD += s.Fracking_EquipmentWear_USDpm3*c_v*(s.Fracking_PressureFactor_ratio+(s.Fracking_PressureFactor_scale*Pstim)**3.0)
    
    #result
    s.Cost_Stimulation_h = Stimulation_Time_h
    s.Cost_Stimulation_USD = -Stimulation_Cost_USD
    return Stimulation_Time_h, -Stimulation_Cost_USD

#**********************************************************************************************************
### calcuate facilities costs
#**********************************************************************************************************
def facilities(s=[],b=[],visuals=False):
    #simultaneous surface facilities with sequential drilling and stimulation
    Construction_Time_yr = s.Cost_ConstructionTime_yr
    
    #use maximum power model (compare to running power model)
    if np.max(b) > 0:
        Construction_Cost_USD = np.max(b)*s.Cost_EquipmentCost_USDpkW*(0.6+3.0*(np.max(b)*1e-3)**(-0.5))
    else:
        Construction_Cost_USD = 0.0
    
    #result
    s.Cost_Facilities_yr = Construction_Time_yr
    s.Cost_Facilities_USD = -Construction_Cost_USD
    return Construction_Time_yr, -Construction_Cost_USD

#**********************************************************************************************************
### total seismic costs (hydrofracs)
#**********************************************************************************************************
def seismic(s=[],visuals=False):
    Mwmax = s.Seismicity_MaxQuake_Mw
    s.Cost_Seismic_USD = s.Cost_Seismic_USDpMw*np.exp(Mwmax*s.Cost_Seismic_exp)
    return s.Cost_Seismic_USD

#**********************************************************************************************************
### total capital costs (excludes planning)
#**********************************************************************************************************
def capital(s=[],visuals=False):
    #simultaneous surface facilities with sequential drilling and stimulation
    Capital_Time_yr = s.Cost_PlanningTime_yr + np.max([s.Cost_Facilities_yr,(s.Cost_Drilling_h+s.Cost_Stimulation_h)/(364.75*24)])
    
    #grants, drilling, stimulation, facilities, and grid connection
    Capital_Cost_USD = s.Cost_Grants_USD
    Capital_Cost_USD += s.Cost_Drilling_USD
    Capital_Cost_USD += s.Cost_Stimulation_USD
    Capital_Cost_USD += s.Cost_Facilities_USD
    Capital_Cost_USD += s.Cost_Gridconnect_USD
    Capital_Cost_USD += s.Cost_Seismic_USD
    
    #result
    s.Cost_Capital_yr = Capital_Time_yr
    s.Cost_Capital_USD = Capital_Cost_USD
    return Capital_Time_yr, Capital_Cost_USD    

#**********************************************************************************************************
### calculate cashflow from timeseries information
#**********************************************************************************************************
# def evaluate(s=iogt.setup(),w=wells.gen_wells(iogt.setup()),b=[1000.0],visuals=False):
def evaluate(s,w,ts,ther,bulk,pump,heat,visuals=False):
    #capital costs
    time, cost = drilling(s,w,visuals)
    time, cost = fracking(s,w,visuals)
    time, cost = facilities(s,bulk,visuals)
    cost = seismic(s)
    time, cost = capital(s)
    
    #match fields from setup to cashflow
    econ = cashflow()
    for item in vars(econ):
        if item in vars(s):
            setattr(econ,item,getattr(s,item))
    
    #power with matched timelines
    gap = int(1.0+s.Cost_Capital_yr)
    projectlife = gap + int(s.Cost_Lifespan_yr)
    econ.Power_Gross_kWh = np.zeros(projectlife)
    econ.Power_Thermal_kWh = np.zeros(projectlife)
    econ.Power_Pumping_kWh = np.zeros(projectlife)
    econ.Power_Heat_kWh = np.zeros(projectlife)
    y = 0
    i0 = 0
    for i in range(0,len(ts)):
        #check if year has advanced (relevant when timesteping is fine)
        if int((ts[i]+1.0)/yr) > y:
            dt = (ts[i]-ts[i0])/(60*60)
            econ.Power_Heat_kWh[gap+y] = 0.5*(ther[i]+ther[i0])*dt
            econ.Power_Gross_kWh[gap+y] = 0.5*(bulk[i]+bulk[i0])*dt
            econ.Power_Pumping_kWh[gap+y] = 0.5*(pump[i]+pump[i0])*dt
            econ.Power_Thermal_kWh[gap+y] = 0.5*(heat[i]+heat[i0])*dt
            y += 1
            i0 = i
            
    #cap thermal production
    for i in range(0,len(econ.Power_Thermal_kWh)):
        if econ.Power_Thermal_kWh[i] > s.Sales_HeatDemand_kWhpyr:
            econ.Power_Thermal_kWh[i] = s.Sales_HeatDemand_kWhpyr
    #compute
    econ.re_init()
    
    #store key results with match fields
    for item in vars(econ):
        if np.isscalar(getattr(econ,item)):
            setattr(s,item,getattr(econ,item))
    
    #visuals
    if visuals:
        yrs = len(econ.Cost_Capex_USD)
        xs = np.linspace(0,yrs,yrs+1)+2025
        fig = pylab.figure(figsize=(10.0,8.0),dpi=100)
        ax = fig.add_subplot(211)
        ax.axhline(y=0.0,linestyle='--',color='grey',linewidth=1)
        ax.plot(xs,[0]+list(econ.Cost_Capex_USD),'-',label='CAPEX',color='blue')
        ax.plot(xs,[0]+list(econ.Profit_BeforeTax_USD),'-',label='PBT',color='green')
        ax.plot(xs,[0]+list(econ.Profit_AfterTax_USD),'-',label='PAT',color='orange')
        ax.plot(xs,[0]+list(econ.Profit_Equity_USD),'-',label='EQUITY',color='magenta')
        ax.legend(loc='upper left',ncol=2,fontsize=10)
        ax.set_ylabel('Cash Flow (USD)',fontsize=10)
        ax.set_xlabel('Date (yr)',fontsize=10)
        ax.yaxis.set_major_formatter(pylab.matplotlib.ticker.StrMethodFormatter('${x:,.0f}'))
        pylab.tight_layout()
        
        ax2 = fig.add_subplot(212)
        ax2.axhline(y=0.0,linestyle='--',color='grey',linewidth=1)
        # IRR10 = (econ.Cost_Capital_USD*s.Equity_Contribution_ratio+econ.Cost_Planning_USD)*(2.0-(1.0+.1-econ.Tax_Inflation_ratio)**np.linspace(0,yrs,yrs+1))
        # IRR10 = list(np.zeros(2+int(econ.Cost_Capital_yr-0.001))) + list(IRR10[0:yrs-1-int(econ.Cost_Capital_yr-0.001)])
        # print(IRR10)
        # print(len(IRR10))
        # ax2.plot(xs,IRR10,':',label='IRR10',color='silver')
        def vc(vals):
            out = np.zeros(yrs+1)
            for i in range(0,yrs):
                out[i+1] = np.sum(vals[0:i])
            return out
        ax2.plot(xs,vc(econ.Cost_Capex_USD),'-',label='CAPEX',color='blue')
        ax2.plot(xs,vc(econ.Profit_BeforeTax_USD),'-',label='PBT',color='green')
        ax2.plot(xs,vc(econ.Profit_AfterTax_USD),'-',label='PAT',color='orange')
        ax2.plot(xs,vc(econ.Profit_Equity_USD),'-',label='EQUITY',color='magenta')
        ax2.legend(loc='upper left',ncol=2,fontsize=10)
        ax2.set_ylabel('Cumulative Value (USD)',fontsize=10)
        ax2.set_xlabel('Date (yr)',fontsize=10)
        ax2.yaxis.set_major_formatter(pylab.matplotlib.ticker.StrMethodFormatter('${x:,.0f}'))
        pylab.tight_layout()
    
    #result
    return econ

#**********************************************************************************************************
### code verification (Texas)
#**********************************************************************************************************
if False:
    #**********************************************************************************************************
    ### costs for triplet well at Durham 6500m depth 3500 m lateral 8.625in diameter surface, intermediate, openhole
    #**********************************************************************************************************
    #setup
    s = iogt.setup()
    s.Strategy_Design_type = 'DCM'
    s.Strategy_Target_type = 'DEEP'
    s.Strategy_Stim_type = 'RADI'
    s.Well_Pattern_count = 1 #pattern
    s.Well_Spacing_m = 150.0 #m
    s.Well_DeviatedLength_m = 1500.0 #m
    s.Well_ProducerProportion_ratio = 1.0 #m/m
    s.Well_Azimuth_rad = 90.0*deg #rad
    s.Well_Dip_rad = 0.0*deg #rad
    s.Well_RotationPhase_rad = 0.0*deg #rad
    s.Well_RotationToe_rad = 0.0*deg #rad
    s.Well_RotationSkew_rad = 0.0*deg #rad
    s.Well_ProductionDiameter_m = 0.0254*7.0 #8.75 #m
    s.Stress_ShminToSv_ratio = 0.8
    s.Stimulation_InjectionRate_m3ps = 0.0927 #m3ps
    s.Stimulation_ProppantConcentration_m3pm3 = 0.059 #m3/m3
    s.Sales_HeatDemand_kWhpyr = 100578.8571 #kWh/yr
    s.Circulation_InjectionTemperature_K = 60.0+273.15 #K
    s.Equity_Contribution_ratio = 0.4
    s.Sales_ContractTime_yr = 30
    s.Debt_Tenor_yr = 30.0
    
    #tuning
    s.Stimulation_Leakoff_ratio = 2.40 #3.55 #m3/m3
    #Drilling Cost Model
    s.Drilling_TimeMobilize_h = 0.0 #336.0 #hr
    s.Drilling_CostMobilize_USD = 860000.0+1600000.0 #1500000.0 #USD
    s.Drilling_TimeRigWalk_h = 48.0 #hr
    s.Drilling_CostRigWalk_USD = 250000.0 #USD
    # s.Drilling_TimeCut_hpm4 = 6.6e-9 #2.3e-11 #h/m4
    # s.Drilling_TimeCut_hpm3 = 0.020 #0.5 #h/m3
    s.Drilling_TimeCut_hpm2 = 6.6e-9 #2.3e-11 #h/m4
    s.Drilling_TimeCut_hpm = 0.020 #0.5 #h/m3
    s.Drilling_CostCut_USDph = 4000.0 #5082.1 #$/hr
    s.Drilling_CostCut_USDpm4 = 1.0e-5 #1.5e-5 #4.5e-5 #$/m4
    s.Drilling_CostCut_USDpm3 = 1500.0 #500.0 #1000.0 #$/m3
    # s.Drilling_TimeCasing_hpm4 = 0.00 #2.0e-12 #h/m4
    # s.Drilling_TimeCasing_hpm3 = 0.01 #0.15 #h/m3
    # s.Drilling_TimeCasing_hpm2 = 0.00 #2.0e-12 #h/m4
    s.Drilling_TimeCasing_hpm = 0.01 #0.15 #h/m3
    s.Drilling_CostCasing_USDph = 3500.0 #4001.3 #$/hr
    s.Drilling_CostCasing_USDpm4 = 0.5 #1.2 #1.4 #$/m4
    s.Drilling_CostCasing_USDpm3 = 1500.0 #2000.0 #1460.0 #$/m3
    s.Drilling_Learning_ratio = 0.25 #0.20 #learning rate
    
    #Stimulation Cost Model
    s.Fracking_Fixed_USD = 500000.0 #1233149.0 #$
    s.Fracking_Depth_USD = 1000.0 #638.0 #$/m
    s.Fracking_Hourly_USDph = 5000.0 #5973.3 #$/hr
    s.Fracking_Sand_USDpkg = 0.662 #$/kg #ISP
    s.Fracking_Fuel_USDpkWh = 0.15 #0.451 #$/kWh
    s.Fracking_Water_USDpm3 = 22.2 #$/m3
    s.Fracking_PressureFactor_ratio = 0.8
    s.Fracking_PressureFactor_scale = 1e-8 #Pa
    s.Fracking_EquipmentWear_USDpm3 = 9.85 #$/m3
    s.Fracking_TimeFactor_hph = 1.95
    
    #zero depth vtk at ground level
    def zero_depth(setup,well):
        for i in well:
            i.c0 += np.asarray([0.0, 0.0, setup.Well_TargetDepth_m])
            i.c1 += np.asarray([0.0, 0.0, setup.Well_TargetDepth_m])
    
    #header
    print('*** Validation reference: GeoDT - cost models 10-14-25.xlsx ***')
    root = 'Colorado'
    case = 0

    #vertical twin (GLADE)
    case += 1
    s.Well_TargetDepth_m = 6100.0-1500.0 #m
    s.Well_Dip_rad = 80*deg #rad
    s.Well_DeviatedLength_m = 3000.0
    s.Well_RotationPhase_rad = 90*deg #rad
    s.Well_Spacing = 30.0 #m
    s.Well_ProducerCount_wells = 1 #wells
    s.Stimulation_TargetClusters_count = 0 #clusters
    s.Stimulation_TargetRadius_m = 150.0 #m
    s.re_init()
    name = '%s_%.0f_%.0f_%i_%.0f_%.0f' %(root,s.Well_TargetDepth_m,s.Well_DeviatedLength_m,s.Well_ProducerCount_wells+s.Well_InjectorCount_wells,
                                        s.Stimulation_TargetRadius_m,s.Stimulation_TargetClusters_count)
    w = wells.gen_wells(s)
    for i in w:
        if i.type in ['perfcluster','procluster']:
            i.type = 'screen'
    wells.export(name+'.csv',w)
    time, cost = drilling(s,w,True)
    time, cost = fracking(s,w,True)
    zero_depth(s,w)
    wells.vtk(w,name)
    pylab.tight_layout()
    pylab.savefig(name+'.png', format='png')
    print('\n Case %i: %s' %(case, name))
    print('    GeoDT: Drill = $%.0f and %.0f days; Stim = $%.0f, %.0f days, and %.0f m3 water' 
          %(s.Cost_Drilling_USD, s.Cost_Drilling_h/24.0, 
            s.Cost_Stimulation_USD,s.Cost_Stimulation_h/24.0,
            s.Stimulation_TargetVolume_m3pfrac*s.Stimulation_TargetClusters_count*s.Well_InjectorCount_wells))
    print('    PASON: Drill = $%.0f and %.0f days; Stim = $%.0f, %.0f days, and %.0f m3 water' 
          %(-18000000, 48.0,
            0.0, 0.0,
            0.0))
    
    #base case
    root = 'Texas'
    case += 1
    s.Well_TargetDepth_m = 6500.0 #m
    s.Well_Dip_rad = 0.0*deg
    s.Well_RotationPhase_rad = 0.0*deg
    s.Well_DeviatedLength_m = 1500.0
    s.Well_ProducerCount_wells = 1 #wells
    s.Stimulation_TargetClusters_count = 27 #clusters
    s.Stimulation_TargetRadius_m = 150.0 #m
    s.re_init()
    name = '%s_%.0f_%.0f_%i_%.0f_%.0f' %(root,s.Well_TargetDepth_m,s.Well_DeviatedLength_m,s.Well_ProducerCount_wells+s.Well_InjectorCount_wells,
                                        s.Stimulation_TargetRadius_m,s.Stimulation_TargetClusters_count)
    w = wells.gen_wells(s)
    wells.export(name+'.csv',w)
    time, cost = drilling(s,w,True)
    time, cost = fracking(s,w,True)
    zero_depth(s,w)
    wells.vtk(w,name)
    pylab.tight_layout()
    pylab.savefig(name+'.png', format='png')
    print('\n Case %i: %s' %(case, name))
    print('    GeoDT: Drill = $%.0f and %.0f days; Stim = $%.0f, %.0f days, and %.0f m3 water' 
          %(s.Cost_Drilling_USD, s.Cost_Drilling_h/24.0, 
            s.Cost_Stimulation_USD,s.Cost_Stimulation_h/24.0,
            s.Stimulation_TargetVolume_m3pfrac*s.Stimulation_TargetClusters_count*s.Well_InjectorCount_wells))
    print('    XLSX: Drill = $%.0f and %.0f days; Stim = $%.0f, %.0f days, and %.0f m3 water' 
          %(-28800000, 163,
            -17567294, 20.0,
            68583.4))
    
    #more wells
    case += 1
    s.Well_TargetDepth_m = 6500.0 #m
    s.Well_DeviatedLength_m = 1500.0
    s.Well_ProducerCount_wells = 5 #wells
    s.Stimulation_TargetClusters_count = 27 #clusters
    s.Stimulation_TargetRadius_m = 150.0 #m
    s.re_init()
    name = '%s_%.0f_%.0f_%i_%.0f_%.0f' %(root,s.Well_TargetDepth_m,s.Well_DeviatedLength_m,s.Well_ProducerCount_wells+s.Well_InjectorCount_wells,
                                        s.Stimulation_TargetRadius_m,s.Stimulation_TargetClusters_count)
    w = wells.gen_wells(s)
    wells.export(name+'.csv',w)
    time, cost = drilling(s,w,True)
    time, cost = fracking(s,w,True)
    zero_depth(s,w)
    wells.vtk(w,name)
    pylab.tight_layout()
    pylab.savefig(name+'.png', format='png')
    print('\n Case %i: %s' %(case, name))
    print('    GeoDT: Drill = $%.0f and %.0f days; Stim = $%.0f, %.0f days, and %.0f m3 water' 
          %(s.Cost_Drilling_USD, s.Cost_Drilling_h/24.0, 
            s.Cost_Stimulation_USD,s.Cost_Stimulation_h/24.0,
            s.Stimulation_TargetVolume_m3pfrac*s.Stimulation_TargetClusters_count*s.Well_InjectorCount_wells))
    print('    XLSX: Drill = $%.0f and %.0f days; Stim = $%.0f, %.0f days, and %.0f m3 water' 
          %(-104100000.0, 505,
            -17567294*4, 20.0*4,
            68583.4*4))
    
    #shallower
    case += 1
    s.Well_TargetDepth_m = 3500.0 #m
    s.Well_DeviatedLength_m = 1500.0
    s.Well_ProducerCount_wells = 1 #wells
    s.Stimulation_TargetClusters_count = 27 #clusters
    s.Stimulation_TargetRadius_m = 150.0 #m
    s.re_init()
    name = '%s_%.0f_%.0f_%i_%.0f_%.0f' %(root,s.Well_TargetDepth_m,s.Well_DeviatedLength_m,s.Well_ProducerCount_wells+s.Well_InjectorCount_wells,
                                        s.Stimulation_TargetRadius_m,s.Stimulation_TargetClusters_count)
    w = wells.gen_wells(s)
    wells.export(name+'.csv',w)
    time, cost = drilling(s,w,True)
    time, cost = fracking(s,w,True)
    zero_depth(s,w)
    wells.vtk(w,name)
    pylab.tight_layout()
    pylab.savefig(name+'.png', format='png')
    print('\n Case %i: %s' %(case, name))
    print('    GeoDT: Drill = $%.0f and %.0f days; Stim = $%.0f, %.0f days, and %.0f m3 water' 
          %(s.Cost_Drilling_USD, s.Cost_Drilling_h/24.0, 
            s.Cost_Stimulation_USD,s.Cost_Stimulation_h/24.0,
            s.Stimulation_TargetVolume_m3pfrac*s.Stimulation_TargetClusters_count*s.Well_InjectorCount_wells))
    print('    XLSX: Drill = $%.0f and %.0f days; Stim = $%.0f, %.0f days, and %.0f m3 water' 
          %(-20853000.0, 115,
            -17567294, 20.0,
            68583.4))
    
    #deeper
    case += 1
    s.Well_TargetDepth_m = 7400.0 #m
    s.Well_DeviatedLength_m = 1500.0
    s.Well_ProducerCount_wells = 1 #wells
    s.Stimulation_TargetClusters_count = 27 #clusters
    s.Stimulation_TargetRadius_m = 150.0 #m
    s.re_init()
    name = '%s_%.0f_%.0f_%i_%.0f_%.0f' %(root,s.Well_TargetDepth_m,s.Well_DeviatedLength_m,s.Well_ProducerCount_wells+s.Well_InjectorCount_wells,
                                        s.Stimulation_TargetRadius_m,s.Stimulation_TargetClusters_count)
    w = wells.gen_wells(s)
    wells.export(name+'.csv',w)
    time, cost = drilling(s,w,True)
    time, cost = fracking(s,w,True)
    wells.vtk(w,name)
    pylab.tight_layout()
    pylab.savefig(name+'.png', format='png')
    print('\n Case %i: %s' %(case, name))
    print('    GeoDT: Drill = $%.0f and %.0f days; Stim = $%.0f, %.0f days, and %.0f m3 water' 
          %(s.Cost_Drilling_USD, s.Cost_Drilling_h/24.0, 
            s.Cost_Stimulation_USD,s.Cost_Stimulation_h/24.0,
            s.Stimulation_TargetVolume_m3pfrac*s.Stimulation_TargetClusters_count*s.Well_InjectorCount_wells))
    print('    XLSX: Drill = $%.0f and %.0f days; Stim = $%.0f, %.0f days, and %.0f m3 water' 
          %(-30900000, 179,
            -17567294, 20.0,
            68583.4))