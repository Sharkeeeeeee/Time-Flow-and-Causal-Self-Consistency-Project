from decimal import Decimal,ROUND_HALF_UP
import json,math
def cost(tasks,shots):return Decimal(tasks)*Decimal('.30')+Decimal(shots)*Decimal('.00145')
def money(x):return str(x.quantize(Decimal('.01'),rounding=ROUND_HALF_UP))
paired=4*math.sqrt(math.log(2*4*2400/.01)/(2*4096))
needed=math.ceil(16*math.log(2*4*2400/.01)/(2*.1**2))
assert paired>.1 and needed==11575
result={'status':'PASS_PLANNING_ARITHMETIC','paired_radius':paired,'target':.1,
        'zero_observed_difference_cannot_pass':paired>.1,'necessary_width_shots':needed,'power_guaranteed':False,
        'prices_usd':{'small_two_epoch_example':money(cost(34,34*4096)),
        'one_full_calibration_example':money(cost(162,162*2048)),
        'six_models_three_independent_calibrations':money(18*cost(162,162*2048))},
        'rates_checked':'2026-10-14','pricing_source':'https://aws.amazon.com/braket/pricing/',
        'hardware_executions':0,'external_review_established':False}
assert result['prices_usd']=={'small_two_epoch_example':'212.13','one_full_calibration_example':'529.68','six_models_three_independent_calibrations':'9534.15'}
print(json.dumps(result,indent=2))
