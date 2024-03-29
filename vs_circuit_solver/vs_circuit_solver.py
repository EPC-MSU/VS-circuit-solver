# vs_circuit_solver.py
# версия 0.1
# язык Python
#
# программа подбора значений R,C для вариантов электронной схемы
# исходя из моделирования подобной схемы в ngspice
# поставляется без всякой оптимизации, ибо имеет целью установление методики
# расчета таких вещей и определения границ применимости этой методики
#
# автор В.Симонов, 22-июль-2020
# vasily_simonov@mail.ru, github.com/vasily84
#
# license : это модуль в любом виде можно использовать в любых целях.
# Ссылка на автора приветствуется, но не является обязательной
#


import scipy.optimize as spo
import scipy.fft as scf
import math
import numpy as np
import matplotlib.pyplot as plt
from ctypes import c_double
import json
# внешние модули
import MySpice.MySpice as spice
import ivcmp.ivcmp as ivcmp
import gc


# SETTINGS ################################################################

# метод сравнения кривых тока и напряжения
# может быть : 'ivcmp','type_ps'
MISFIT_METHOD = 'ivcmp'
# MISFIT_METHOD = 'type_ps'


# частота, Гц
INIT_F = 1e4
# амплитудное напряжение, Вольт, может изменится при загрузке внешнего файла данных
INIT_V = 5

# токоограничивающий резистор, Ом
INIT_Rcs = 1e2

# SIGNAL/NOISE ratio
INIT_SNR = 120.0
# INIT_SNR = 35.0

# число циклов колебаний напряжения в записи
INIT_CYCLE = 10

# падение напряжения на диоде
# Диод считается полностью проводимым при напряжении больше чем DIODE_VOLTAGE,
# при меньшем полность закрыт. (Приближение)
DIODE_VOLTAGE = 0.7

#
SMALL_VOLTAGE = 0.1

# "огромное сопротивление".
HUGE_R = 1e10  # 10 ГОм

# "большое сопротивление"
BIG_R = 1e8  # 100 МОм

# "мизерное сопротивление"
NULL_R = 1e-6  # 1 мкОм

# "мизерная емкость","огромная емкость"
NONE_C = 1e-15  # 0.001 пФ
HUGE_C = 1e-3  # 1000 мкФ
# погрешность подбора кривых- критерий остановки. Подбор длится до тех пор,
# пока функция сравнения не вернет значение CompareIvc()<=IVCMP_TOLERANCE
# IVCMP_TOLERANCE = 5e-3
IVCMP_TOLERANCE = 6e-2


# погрешность подбора номиналов в процентах. Номиналы емкостей считаются по
# реактивному сопротивлению!. Подробности см. scipy.minimize(method='Powell')
VALUES_TOLERANCE = 1e-2

# число вычислений функции в процессе оптимизации. При малых значениях-
# минимально возможное число
MAXFEV = 100

# число точек в массивах тока и напряжения, может измениться при загрузке
# внешнего файла данных
MAX_NUM_POINTS = 100


min_ivc = 1
#############################################################################

# результат последнего моделирования в PySpice
analysis = None

# целевая кривая с током. Та, которую мы подбираем
target_VCurrent = None
# измеренное прибором напряжение в точке после резистора Rcs
target_input_dummy = None

target_fileName = ''

# целевая кривая с током для сравнения в библиотеке ivcmp
target_IVCurve = None


# название временного файла схемы для запуска PySpice
circuit_SessionFileName = 'var1.cir'

# список значений для файла шаблона схемы. Число элементов - не меньше, чем
# знаков {} в файле шаблона схемы
# Xi_long = [0.,0.,0.,0., 0.,0.,0., 0.,0.,0.,0.]
Xi_long = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])


# Маска оптимизируемых параметров - список булевого типа, например -
# Xi_long = [a, b, c, d]
# Xi_mask = [False,True,False,True] -> X_short = [b,d]
Xi_mask = [False, False, False, False, False, False, False, False, False, False, False]


# ФУНКЦИИ ДЛЯ ШАБЛОНА, ЦЕЛЕВОЙ МОДЕЛИ И МАСКИ ПАРАМЕТРОВ ##################
def Xi_unroll(x_short):
    XL = Xi_long.copy()

    j = 0
    for i in range(0, len(Xi_mask)):
        if Xi_mask[i]:
            XL[i] += x_short[j]
            j += 1

    return XL


def Xi_pack(Xi_):
    xi = []

    for i in range(0, len(Xi_mask)):
        if Xi_mask[i]:
            xi += [Xi_[i]]

    return xi


# установить все известные номиналы
def set_circuit_nominals(nominals):
    global Xi_long
    Xi_long = nominals.copy()


def reset_Xi_variable():
    for i in range(len(Xi_mask)):
        Xi_mask[i] = False


def set_Xi_variable(vlist):
    for v in vlist:
        if v == 'R1':
            Xi_mask[0] = True
        if v == 'C1':
            Xi_mask[1] = True
        if v == '_R_C1':
            Xi_mask[2] = True
        if v == '_R_D1':
            Xi_mask[3] = True

        if v == 'R2':
            Xi_mask[4] = True
        if v == 'C2':
            Xi_mask[5] = True
        if v == '_R_C2':
            Xi_mask[6] = True

        if v == 'R3':
            Xi_mask[7] = True
        if v == 'C3':
            Xi_mask[8] = True
        if v == '_R_C3':
            Xi_mask[9] = True
        if v == '_R_D3':
            Xi_mask[10] = True


def sign(value):
    if value < 0:
        return -1
    else:
        return 1


# инициализировать целевую модель, промоделировав файл схемы
def init_target_by_circuitFile(fileName=circuit_SessionFileName):
    global target_VCurrent, target_input_dummy, target_IVCurve
    global circuit_SessionFileName
    global Z123_sch
    Z123_sch = None

    var1 = circuit_SessionFileName
    circuit_SessionFileName = fileName
    process_circuitFile()
    circuit_SessionFileName = var1

    target_VCurrent = analysis.VCurrent
    target_input_dummy = analysis.input_dummy

    iv_curve = ivcmp.IvCurve()
    iv_curve.length = MAX_NUM_POINTS-1
    for i in range(MAX_NUM_POINTS-1):
        iv_curve.voltages[i] = c_double(analysis.input_dummy[i])  # Ток и напряжение были поменяны местами
        iv_curve.currents[i] = c_double(analysis.VCurrent[i])

    min_var_c = 0.01 * np.max(iv_curve.currents[:MAX_NUM_POINTS-1])  # value of noise for current
    min_var_v = 0.01 * np.max(iv_curve.voltages[:MAX_NUM_POINTS-1])  # value of noise for voltage
    # if (abs(min_var_c) < INIT_V/INIT_Rcs*0.03):
    #    min_var_c = sign(min_var_c)*INIT_V/INIT_Rcs*0.03
    ivcmp.SetMinVC(min_var_v, min_var_c)  # Правильные значения фильтров для корректной работы

    target_IVCurve = iv_curve


# инициализировать целевую модель данными из json файла, установить число точек на кривой MAX_NUM_POINTS
# определенными из файла
def init_target_from_jsnFile(fileName, N):

    global target_fileName

    target_fileName = fileName

    ivc_real = open_board(fileName)
    if ivc_real is None:
        print('open_board() failed')
        return

    print('record number = '+str(N))
    target_voltages = ivc_real["elements"][0]["pins"][N]["iv_curves"][0]["voltages"]
    target_currents = ivc_real["elements"][0]["pins"][N]["iv_curves"][0]["currents"]

    # частота, Гц
    initF = ivc_real["elements"][0]["pins"][N]["iv_curves"][0]["measurement_settings"]["probe_signal_frequency"]
    print('INIT_F = '+str(initF))
    # амплитудное напряжение, Вольт, может изменится при загрузке внешнего файла данных
    initV = ivc_real["elements"][0]["pins"][N]["iv_curves"][0]["measurement_settings"]["max_voltage"]
    print('INIT_V = '+str(initV))
    # токоограничивающий резистор, Ом
    initRcs = ivc_real["elements"][0]["pins"][N]["iv_curves"][0]["measurement_settings"]["internal_resistance"]
    print('INIT_Rcs = '+str(initRcs))
    return initF, initV, initRcs, target_voltages, target_currents


def init_target_Data(target_voltages,
                     target_currents,
                     initF=1e4,
                     initV=5,
                     initRcs=1e2,
                     initSNR=120,
                     cycle=10,
                     ivcmpTolerance=6e-2):
    global MAX_NUM_POINTS, INIT_V, INIT_F, INIT_Rcs
    global target_VCurrent, target_input_dummy, target_IVCurve

    target_input_dummy = target_voltages
    target_VCurrent = target_currents
    INIT_F = initF
    INIT_V = initV
    INIT_Rcs = initRcs
    MAX_NUM_POINTS = len(target_input_dummy)
    print('MAX_NUM_POINTS = '+str(MAX_NUM_POINTS))

    iv_curve1 = ivcmp.IvCurve()
    iv_curve1.length = MAX_NUM_POINTS-1
    for i in range(MAX_NUM_POINTS-1):
        iv_curve1.voltages[i] = c_double(target_input_dummy[i])  # Ток и напряжение были поменяны местами
        iv_curve1.currents[i] = c_double(target_VCurrent[i])

    min_var_c = 0.01 * np.max(iv_curve1.currents[:MAX_NUM_POINTS-1])  # value of noise for current
    min_var_v = 0.01 * np.max(iv_curve1.voltages[:MAX_NUM_POINTS-1])  # value of noise for voltage

    ivcmp.SetMinVC(min_var_v, min_var_c)  # Правильные значения фильтров для корректной работы

    target_IVCurve = iv_curve1
    return


def Xi_to_RC(Xi):
    RC = Xi.copy()

    RC[0] = np.abs(Xi[0])
    RC[1] = np.abs(R_to_C(Xi[1]))  # C1
    RC[2] = np.abs(Xi[2])
    RC[3] = np.abs(Xi[3])

    RC[4] = np.abs(Xi[4])
    RC[5] = np.abs(R_to_C(Xi[5]))  # C2
    RC[6] = np.abs(Xi[6])

    RC[7] = np.abs(Xi[7])
    RC[8] = np.abs(R_to_C(Xi[8]))  # C3
    RC[9] = np.abs(Xi[9])
    RC[10] = np.abs(Xi[10])
    return RC


# в наборе строк шаблона схемы сделать замену {} на значения
# варьирования Xi_values, сохранить заданным с именем
def generate_circuitFile_by_values(Xi_values):
    rc_values = Xi_to_RC(Xi_values)
    with open(circuit_SessionFileName, 'w') as newF:
        newF.write('* cir file corresponding to the equivalent circuit.\n')
        # * Цепь 1
        if rc_values[0] < BIG_R:  # цепь R1 присутствует
            if rc_values[2] >= BIG_R:  # C1 присутствует
                newF.write('R1 _net1 input {:e}\n'.format(rc_values[0]))
                newF.write('C1 _net0 _net1 {:e}\n'.format(rc_values[1]))
            else:  # С1 нет
                newF.write('R1 _net0 input {:e}\n'.format(rc_values[0]))

            if rc_values[3] >= BIG_R:  # D1 присутствует
                newF.write('D1 _net0 0 DMOD_D1 AREA=1.0 Temp=26.85\n')
            else:  # вместо D1 перемычка
                newF.write('R_D1 0 _net0 {:e}\n'.format(rc_values[3]))

        # * Цепь 2
        if rc_values[4] < BIG_R:
            if rc_values[6] >= BIG_R:  # C2 присутствует
                newF.write('R2 _net4 input {:e}\n'.format(rc_values[4]))
                newF.write('C2 0 _net4 {:e}\n'.format(rc_values[5]))
            else:  # вместо С2 перемычка, R2 сразу на землю
                newF.write('R2 0 input {:e}\n'.format(rc_values[4]))

        # * Цепь 3
        if rc_values[7] < BIG_R:
            if rc_values[9] >= BIG_R:  # C3 присутствует
                newF.write('R3 _net3 input {:e}\n'.format(rc_values[7]))
                newF.write('C3 _net2 _net3 {:e}\n'.format(rc_values[8]))
            else:  # С3 нет
                newF.write('R3 _net2 input {:e}\n'.format(rc_values[7]))

            if rc_values[10] >= BIG_R:  # D3 присутствует
                newF.write('D3 0 _net2 DMOD_D1 AREA=1.0 Temp=26.85\n')
            else:  # вместо D3 перемычка
                newF.write('R_D3 0 _net2 {:e}\n'.format(rc_values[10]))

        # есть диоды, добавляем модель
        if (rc_values[10] >= BIG_R) or (rc_values[3] >= BIG_R):
            newF.write('.MODEL DMOD_D1 D (Is=2.22e-10 N=1.65 Cj0=4e-12 M=0.333 '
                       'Vj=0.7 Fc=0.5 Rs=0.0686 Tt=5.76e-09 Ikf=0 Kf=0 Af=1 Bv=75 '
                       'Ibv=1e-06 Xti=3 Eg=1.11 Tcv=0 Trs=0 Ttt1=0 Ttt2=0 Tm1=0 Tm2=0 Tnom=26.85 )\n')

        newF.write('.END')
    # end of with


input_data = None


# промоделировать файл схемы
def process_circuitFile():
    global analysis, input_data

    if input_data is None:
        input_data = spice.Init_Data(INIT_F, INIT_V, INIT_Rcs, INIT_SNR)

    try:
        circuit = spice.LoadFile(circuit_SessionFileName)
    except Exception:
        print('spice.LoadFile() failed.')

    try:
        analysis = spice.CreateCVC1(circuit, input_data, MAX_NUM_POINTS, "input", INIT_CYCLE)
    except Exception:
        print('spice.CreateCVC1() failed.')


# последний анализ перевести в форму, пригодную для сравнения в ivcmp
iv_curve = None


def analysis_to_IVCurve():
    global iv_curve

    if iv_curve is None:
        iv_curve = ivcmp.IvCurve()
        iv_curve.length = MAX_NUM_POINTS-1
    for i in range(MAX_NUM_POINTS-1):
        iv_curve.voltages[i] = c_double(analysis.input_dummy[i])
        iv_curve.currents[i] = c_double(analysis.VCurrent[i])

    return iv_curve


def V_div_I(v, i):
    try:
        r = v/i
    except ArithmeticError:
        r = HUGE_R
    return r


# вывести на график результат моделирования
def analysis_plot(title='', pngName=''):
    plt.figure(1, (20, 10))
    plt.grid()

    # целевая ВАХ
    plt.plot(target_input_dummy, target_VCurrent, color='red')
    # ВАХ результат подбора
    plt.plot(analysis.input_dummy, analysis.VCurrent, color='blue')

    s = ''
    if (not title == ''):
        s = title
    elif not target_fileName == '':
        s = target_fileName

    s = s+', misfit=' + format(misfit_result, '0.5E')+', ivcmp='+format(ivcmp_result, '0.5E')

    plt.title(s)

    plt.xlabel('Напряжение [В]')
    plt.ylabel('Сила тока [А]')

    if (not pngName == ''):
        plt.savefig(pngName)
    plt.xlim([-INIT_V, INIT_V])
    plt.ylim([-INIT_V/INIT_Rcs, INIT_V/INIT_Rcs])
    plt.show()


# ФУНКЦИИ СРАВНЕНИЯ ВАХ ###################################################
def C_to_R(c):
    r = 1/(2.*np.pi*INIT_F*c)
    return r


def R_to_C(r):
    c = 1/(2.*np.pi*INIT_F*r)
    if math.isinf(c):
        c = 1e20
    return c


def analysis_misfit_ivcmp():
    global min_ivc
    step_IVCurve = analysis_to_IVCurve()
    res = ivcmp.CompareIvc(target_IVCurve, step_IVCurve)
    if min_ivc > res:
        min_ivc = res
    return res


# вычислить несовпадение последнего анализа и целевой функции.
def analysis_misfit():
    curr_t = target_VCurrent
    curr_a = analysis.VCurrent
    volt_t = target_input_dummy
    volt_a = analysis.input_dummy

    # метод сравнения кривых по несовпадению кривых мощности.
    # учитывает возможное несогласование фаз сигналов
    if MISFIT_METHOD == 'type_ps':
        fullV_target = np.zeros_like(target_input_dummy)
        fullV_A = np.zeros_like(target_input_dummy)
        signal_target = np.zeros_like(target_input_dummy)
        signal_A = np.zeros_like(target_input_dummy)
        signal_cmp = np.zeros_like(target_input_dummy)

        for i in range(len(fullV_target)):
            # полные напряжения возбуждения
            fullV_target[i] = target_input_dummy[i]+INIT_Rcs*target_VCurrent[i]
            fullV_A[i] = analysis.input_dummy[i]+INIT_Rcs*analysis.VCurrent[i]
            # мощности, ушедшие в нагрузку
            # signal_target[i] = fullV_target[i]*target_VCurrent[i]
            # signal_A[i] = fullV_A[i]*analysis.VCurrent[i]
            signal_target[i] = target_VCurrent[i]
            signal_A[i] = analysis.VCurrent[i]

        # выравнивание фаз по максимуму сигнала
        index_target = np.argmax(fullV_target)
        index_A = np.argmax(fullV_A)
        # фазовый сдвиг в отсчетах
        phase_shift = index_A - index_target + len(signal_target)

        for i in range(len(signal_target)):
            i_A = (i+phase_shift) % len(signal_target)
            # разница мгновенной мощности
            signal_cmp[i] = (signal_target[i]-signal_A[i_A])**2

        return math.fsum(signal_cmp)

    if MISFIT_METHOD == 'power_fft':
        r = scf.rfft(curr_t*volt_t-curr_a*volt_a)
        return math.fsum(r)

    if MISFIT_METHOD == 'sko':
        r = (curr_t-curr_a)
        r2 = np.abs(r)
        return math.fsum(r2)

    if MISFIT_METHOD == 'ivcmp':
        step_IVCurve = analysis_to_IVCurve()
        res = ivcmp.CompareIvc(target_IVCurve, step_IVCurve)
        return res

    ###
    s = "unknown MISFIT_METHOD = '"+str(MISFIT_METHOD)+"'"
    raise RuntimeError(s)


# ФУНКЦИИ РЕШАТЕЛЯ ########################################################
# вектор, обеспечивающий минимум оптимизируемой функции
Xi_result = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])

# текущий найденный минимум оптимизируемой функции
misfit_result = 0.
# результат сравнения найденного минимума по функцией CompareIvc()
ivcmp_result = 0.

# счетчик числа вызовов функции оптимизатором
FitterCount = 0
BestMisfitCount = 0
FITTER_SUCCESS = False


def calculate_misfit(Xi):
    generate_circuitFile_by_values(Xi)
    process_circuitFile()
    misfit = analysis_misfit()
    return misfit


# функция вызывается оптимизатором
def fitter_subroutine(Xargs):
    global Xi_result, misfit_result, FitterCount, BestMisfitCount, ivcmp_result, FITTER_SUCCESS
    FitterCount += 1
    xi = Xi_unroll(Xargs)
    misfit = calculate_misfit(xi)
    if MISFIT_METHOD == 'ivcmp':
        ivcmp_result = misfit
    # первый запуск
    if FitterCount <= 1:
        Xi_result = xi.copy()
        misfit_result = misfit
        ivcmp_result = analysis_misfit_ivcmp()
        BestMisfitCount = 0

    # лучший случай
    if misfit < misfit_result:
        Xi_result = xi.copy()
        misfit_result = misfit
        ivcmp_result = analysis_misfit_ivcmp()
        BestMisfitCount += 1

    # дополнительная проверка
    if ivcmp_result <= IVCMP_TOLERANCE:  # достигли необходимой точности
        FITTER_SUCCESS = True

    return misfit


def fitter_callback(Xk):
    global FITTER_SUCCESS
    if ivcmp_result <= IVCMP_TOLERANCE:  # достигли необходимой точности
        FITTER_SUCCESS = True
        return True

    return False


# запустить автоподбор - сравнение по сумме отклонений точек
def run_fitter(result_cir_file_name='', result_csv_file_name=''):
    Xargs = Xi_pack(Xi_long)

    for i in range(0, len(Xargs)):
        Xargs[i] = 0.

    resX = spo.minimize(fitter_subroutine,
                        Xargs,
                        method='Powell',
                        callback=fitter_callback,
                        options={'maxfev': MAXFEV, 'xtol': VALUES_TOLERANCE})

    if (not result_csv_file_name == ''):
        spice.SaveFile(analysis, result_csv_file_name)
    if (not result_cir_file_name == ''):
        generate_circuitFile_by_values(resX.x)

    return True


# элементарная схема ###
##############################################################################
def Sch_init():
    R = HUGE_R  # 100.
    C = NONE_C  # 1e-6
    sch = {}
    sch['R1'] = R
    sch['C1'] = C  # [1]
    sch['_R_C1'] = NULL_R
    sch['_R_D1'] = HUGE_R

    sch['R2'] = R
    sch['C2'] = C  # [5]
    sch['_R_C2'] = NULL_R

    sch['R3'] = R
    sch['C3'] = C  # [8]
    sch['_R_C3'] = NULL_R
    sch['_R_D3'] = HUGE_R
    return sch


def Sch_get_Xi(sch):
    xi = []
    for k in sch:
        if (k == 'C1') or (k == 'C2') or (k == 'C3'):
            xi += [C_to_R(sch[k])]
        else:
            xi += [sch[k]]

    return xi


def Sch_load_from_Xi(sch, Xi):
    j = 0
    for k in sch:
        if (k == 'C1') or (k == 'C2') or (k == 'C3'):
            sch[k] = R_to_C(Xi[j])
        else:
            sch[k] = Xi[j]
        j += 1


CODE2_COUNT = 4


# ses - сессия варьирования, которую необходимо проинициализировать
# swcode - числовой код,от 0 до 255 включительно, задает положения переключателей
# code2 - дополнительный код, для каждого варианта swcode передавать code2=0,1,2, ...
# до тех пор, пока функция не вернет False
def Session_init_by_approximation(ses, swcode, code2, title=''):
    sch = ses['start_sch']
    res = Z123_approximation(sch, swcode, code2, title)
    Session_set_switchers(ses, swcode)
    Session_run1(ses)

    return res  # функция больше не вызывается


#############################################################################
# ФУНКЦИИ НУЛЕВОГО ПОДБОРА (ПРИСТРЕЛКА) ####################################
#############################################################################

# полное напряжение цепи - до резистора Rcs.
target_fullVoltage = None
# ток цепи, с коррекцией смещения нуля
corrected_VCurrent = None


# ток через нашу упрощенную цепь
def I_from_VR1R2R3(V, R1, R2, R3):
    I_ = V/(INIT_Rcs+R2)
    R1 = np.abs(R1)
    R2 = np.abs(R2)
    R3 = np.abs(R3)
    V2 = R2*I_

    # диод VD1 открыт
    if V2 >= DIODE_VOLTAGE:
        up_part = V*(R1+R2)-R2*DIODE_VOLTAGE
        down_part = R1*R2+R1*INIT_Rcs+R2*INIT_Rcs
        Id = up_part/down_part
        return Id

    # диод VD3 открыт
    if V2 <= -DIODE_VOLTAGE:
        up_part = V*(R3+R2)+R2*DIODE_VOLTAGE
        down_part = R3*R2+R3*INIT_Rcs+R2*INIT_Rcs
        Id = up_part/down_part
        return Id

    # случай, когда диоды VD1 и VD3 закрыты - просто закон ома
    return I_


# сопротивление из известных значений
def R1_from_R2VI(R2, V, I_):
    I2 = V/(INIT_Rcs+R2)
    # диод открыт
    if (I2*R2) < (DIODE_VOLTAGE-SMALL_VOLTAGE):
        print('R1_from_R2VI() Error: DIODE VD1 CLOSED!!')
        raise RuntimeError("R1_from_R2VI() Error: DIODE VD1 CLOSED!!") from None

    up_part = R2*(V - I_*INIT_Rcs - DIODE_VOLTAGE)
    down_part = I_*(R2+INIT_Rcs)-V
    return up_part/down_part


# сопротивление из известных значений
def R3_from_R2VI(R2, V, I_):
    I2 = V/(INIT_Rcs+R2)
    # диод открыт
    if (I2*R2) > -(DIODE_VOLTAGE-SMALL_VOLTAGE):
        raise RuntimeError("R3_from_R2VI() Error: DIODE VD3 CLOSED!!") from None

    up_part = R2*(V-I_*INIT_Rcs+DIODE_VOLTAGE)
    down_part = I_*(R2+INIT_Rcs)-V
    return up_part/down_part


# измерить непосредственно r2
def measure_r2():
    v_r2 = target_input_dummy

    r_summ = 0.
    r_count = 0

    for i in range(len(v_r2)):
        if (np.abs(v_r2[i]) > SMALL_VOLTAGE) and (np.abs(v_r2[i]) < DIODE_VOLTAGE):
            r_i = V_div_I(v_r2[i], corrected_VCurrent[i])
            if (r_i >= HUGE_R):
                continue
            r_summ += np.abs(r_i)
            r_count += 1
    try:
        R = r_summ/r_count
    except ZeroDivisionError:
        R = HUGE_R
    return R


# измерить непосредственно r1
def measure_r1_by_R2(R2):
    i = np.argmax(target_fullVoltage)
    try:
        r = R1_from_R2VI(R2, target_fullVoltage[i], corrected_VCurrent[i])
    except Exception:
        r = NULL_R

    return r


# измерить непосредственно r3
def measure_r3_by_R2(R2):
    i = np.argmin(target_fullVoltage)
    try:
        r = R3_from_R2VI(R2, target_fullVoltage[i], corrected_VCurrent[i])
    except Exception:
        r = NULL_R

    return r


def get_r_high():
    i = np.argmax(target_fullVoltage)
    r = V_div_I(target_fullVoltage[i], corrected_VCurrent[i])
    return r


def get_r_low():
    i = np.argmin(target_fullVoltage)
    r = V_div_I(target_fullVoltage[i], corrected_VCurrent[i])
    return r


def get_r_hight_sub_diode():
    i = np.argmax(target_fullVoltage)
    r = V_div_I(target_fullVoltage[i]-DIODE_VOLTAGE, corrected_VCurrent[i])
    return r


def get_r_low_sub_diode():
    i = np.argmin(target_fullVoltage)
    r = V_div_I(target_fullVoltage[i]+DIODE_VOLTAGE, corrected_VCurrent[i])
    return r


# при вычислений запоминает лучший результат по совпадению
# кривых.
min_r123_misfit = None
min_r123_x = None


def min_r123_subroutine(x):
    global min_r123_misfit, min_r123_x

    r1 = x[0]
    r2 = x[1]
    r3 = x[2]

    E_r123 = np.zeros_like(target_fullVoltage)

    for i in range(len(E_r123)):
        I_ = I_from_VR1R2R3(target_fullVoltage[i], r1, r2, r3)
        E_r123[i] = (I_-corrected_VCurrent[i])**2

    Result = math.fsum(E_r123)
    if (min_r123_misfit is None) or (Result < min_r123_misfit):
        min_r123_x = [r1, r2, r3]
        min_r123_misfit = Result

    return Result


# измерить смещение нуля в пределах напряжений, где диоды закрыты
def measure_zero_drift():
    z_value = 0.
    z_count = 0
    for i in range(len(target_VCurrent)):
        if (np.abs(target_input_dummy[i]) < (DIODE_VOLTAGE)):
            z_value += target_VCurrent[i]
            z_count += 1
    try:
        z_drift = z_value/z_count
    except ZeroDivisionError:
        z_drift = 0.
    print('z_drift='+str(z_drift))
    return z_drift


# проверить границы номиналов емкости,
# установить граничные значения, если выходит за пределы
def C_to_norm(C):
    if C < NONE_C:
        return NONE_C
    if C > HUGE_C:
        return HUGE_C
    return C


def phase_to_norm(phase):
    pass


# инициализация, исходя из того Rcs может быть
# десятки килоОм
Z123_sch = None


def Z123_approximation(sch, swcode, code2, title=''):
    global Z123_sch, target_fullVoltage, min_r123_misfit, corrected_VCurrent

    if Z123_sch is None:
        Z123_sch = Sch_init()

        target_fullVoltage = np.copy(target_VCurrent)
        corrected_VCurrent = np.copy(target_VCurrent)
        for i in range(len(target_input_dummy)):
            target_fullVoltage[i] = target_input_dummy[i]+INIT_Rcs*target_VCurrent[i]
            corrected_VCurrent[i] = target_VCurrent[i]
    else:
        # копирование
        sch['R1'] = Z123_sch['R1']
        sch['C1'] = Z123_sch['C1']
        sch['R2'] = Z123_sch['R2']
        sch['C2'] = Z123_sch['C2']
        sch['R3'] = Z123_sch['R3']
        sch['C3'] = Z123_sch['C3']

        return False  # больше не вызывать

    ########################################################

    r2 = measure_r2()
    r1 = measure_r1_by_R2(r2)
    r3 = measure_r3_by_R2(r2)

    #  обнуляем пристрелку
    min_r123_misfit = None
    # варианты значений сопротивлений схем
    # основной вариант аналитического приближения, срабатывает почти всегда
    min_r123_subroutine([r1, r2, r3])
    # разные варианты с меньшей абсолютной  погрешностью
    # для аналитического приближения
    min_r123_subroutine([r1, HUGE_R, r3])
    min_r123_subroutine([r1, HUGE_R, measure_r3_by_R2(HUGE_R)])
    min_r123_subroutine([measure_r1_by_R2(HUGE_R), HUGE_R, r3])
    min_r123_subroutine([measure_r1_by_R2(HUGE_R), HUGE_R, measure_r3_by_R2(HUGE_R)])
    min_r123_subroutine([measure_r1_by_R2(HUGE_R), HUGE_R, r3])
    min_r123_subroutine([measure_r1_by_R2(HUGE_R), HUGE_R, HUGE_R])
    min_r123_subroutine([HUGE_R, HUGE_R, measure_r3_by_R2(HUGE_R)])
    # разные варианты с меньшей абсолютной  погрешностью
    # для приближения диода с идеальной ВАХ
    r1_0 = get_r_high()
    r1_d = get_r_hight_sub_diode()
    r3_0 = get_r_low()
    r3_d = get_r_low_sub_diode()

    min_r123_subroutine([r1_d, HUGE_R, r3_d])
    min_r123_subroutine([r1_d, r3_0, HUGE_R])
    min_r123_subroutine([HUGE_R, r1_0, r3_d])
    min_r123_subroutine([HUGE_R, r3_0, HUGE_R])
    min_r123_subroutine([HUGE_R, r1_0, HUGE_R])
    min_r123_subroutine([NULL_R, r2, NULL_R])
    min_r123_subroutine([NULL_R, r2, r3])
    min_r123_subroutine([r1, r2, NULL_R])
    # маловероятно, но пусть будет
    min_r123_subroutine([r1, NULL_R, r3])

    r1 = np.abs(min_r123_x[0])
    r2 = np.abs(min_r123_x[1])
    r3 = np.abs(min_r123_x[2])

    Rc1 = 1./(1./r1+1./r2)
    Rc2 = 1./(1./r1+1./r2+1./r3)
    Rc3 = 1./(1./r2+1./r3)

    phase_1 = 360*(np.argmax(target_fullVoltage)-np.argmax(target_VCurrent))/MAX_NUM_POINTS
    phase_3 = 360*(np.argmin(target_fullVoltage)-np.argmin(target_VCurrent))/MAX_NUM_POINTS
    print('phase_1='+str(phase_1))
    print('phase_3='+str(phase_3))
    phase_1 = np.abs(phase_1) % 90
    phase_3 = np.abs(phase_3) % 90

    if phase_1 < 5:
        phase_1 = 5
    if phase_3 < 5:
        phase_3 = 5
    if phase_1 > 85:
        phase_1 = 85
    if phase_3 > 85:
        phase_3 = 85

    phase_2 = (phase_1+phase_3)/2.
    print('phase_1*='+str(phase_1))
    print('phase_2*='+str(phase_2))
    print('phase_3*='+str(phase_3))

    с1 = R_to_C(Rc1*np.cos(phase_1*np.pi/180))
    с1 = C_to_norm(с1)
    с2 = R_to_C(Rc2*np.cos(phase_2*np.pi/180))
    с2 = C_to_norm(с2)
    с3 = R_to_C(Rc3*np.cos(phase_3*np.pi/180))
    с3 = C_to_norm(с3)

    Z123_sch['R1'] = r1
    Z123_sch['C1'] = с1
    Z123_sch['R2'] = r2
    Z123_sch['C2'] = с2
    Z123_sch['R3'] = r3
    Z123_sch['C3'] = с3

    str_0 = '\nr1_o={:2.1e}, r2_o={:2.1e}, r3_o={:2.1e}'.format(r1, r2, r3)
    plt.title('Пристрелка '+title+str_0)
    plt.plot(target_input_dummy, target_VCurrent, c='red')

    print('r1_o = '+str(r1))
    print('r2_o = '+str(r2))
    print('r3_o = '+str(r3))
    print('с1_o = '+str(с1))
    print('с2_o = '+str(с2))
    print('с3_o = '+str(с3))

    curr_r123 = np.zeros_like(target_fullVoltage)
    for i in range(len(curr_r123)):
        curr_r123[i] = I_from_VR1R2R3(target_fullVoltage[i], r1, r2, r3)

    plt.plot(target_input_dummy, curr_r123, c='blue')
    plt.legend(['реальные даные', 'Н.У. подбора'])
    plt.show()

    # plt.plot(target_input_dummy)
    # plt.show()
    # plt.plot(target_VCurrent)
    # plt.show()

    # именно такое копирование, ибо надо сохранить ссылку
    sch['R1'] = Z123_sch['R1']
    sch['C1'] = Z123_sch['C1']
    sch['R2'] = Z123_sch['R2']
    sch['C2'] = Z123_sch['C2']
    sch['R3'] = Z123_sch['R3']
    sch['C3'] = Z123_sch['C3']

    return

#############################################################################
#############################################################################
#############################################################################


def Sch_saveToFile(sch, fileName):
    global circuit_SessionFileName
    s = circuit_SessionFileName
    circuit_SessionFileName = fileName
    try:
        Session_run1(sch)
    except Exception:
        with open(fileName, 'w') as newF:
            json.dump(sch, newF)

    print(sch['Xi_variable'])
    circuit_SessionFileName = s
    return


def init_target_by_Sch(sch):
    global Z123_sch
    Z123_sch = None

    generate_circuitFile_by_values(Sch_get_Xi(sch))
    init_target_by_circuitFile()
    return


#############################################################################
def Session_create(start_sch):
    s = {}
    s['start_sch'] = start_sch
    return s


# выполнить схему один раз
def Session_run1(session):
    global misfit_result, ivcmp_result
    try:
        sch = session['result_sch']
    except KeyError:
        sch = session['start_sch']

    xi = Sch_get_Xi(sch)
    set_circuit_nominals(xi)
    session['misfit'] = calculate_misfit(xi)
    misfit_result = session['misfit']
    if MISFIT_METHOD == 'ivcmp':
        ivcmp_result = misfit_result


# запустить подбор для сессии
def Session_run_fitter(session):
    global FitterCount
    FitterCount = 0
    try:
        sch = session['result_sch']
    except KeyError:
        sch = session['start_sch']
    else:
        session['start_sch'] = sch

    var_list = session['Xi_variable']
    set_circuit_nominals(Sch_get_Xi(sch))
    set_Xi_variable(var_list)

    try:
        run_fitter()
    except Exception:
        print('NGSPICE EXCEPTION')

    sch2 = Sch_init()
    Sch_load_from_Xi(sch2, Xi_result)
    session['result_sch'] = sch2

    session['misfit'] = misfit_result
    session['fCount'] = FitterCount
    session['mCount'] = BestMisfitCount


# проверить, имеет ли смысл такая установка переключателей в схеме
def is_valid_switchers(swcode):
    if swcode & (1+2+3):  # все ветви заглушены
        return False

    # все разыгрывание по заглушенной первой ветке
    if (swcode == 1) or (swcode == 1+8) or (swcode == 1+16) or (swcode == 1+8+16):
        return False

    # все разыгрывание по заглушенной второй ветке
    if (swcode == 2) or (swcode == 2+128):
        return False

    # все разыгрывание по заглушенной третьей ветке
    if (swcode == 3) or (swcode == 3+32) or (swcode == 3+64) or (swcode == 3+32+64):
        return False

    return True


# установить переключатели для схемы.
def Session_set_switchers(session, swcode):
    sch = session['start_sch']
    var_list = []

    if swcode & 1:  # ветка 1
        sch['R1'] = HUGE_R
    else:
        var_list += ['R1']

    if swcode & 2:  # ветка 2
        sch['R2'] = HUGE_R
    else:
        var_list += ['R2']

    if swcode & 4:  # ветка 3
        sch['R3'] = HUGE_R
    else:
        var_list += ['R3']

    if swcode & 8:  # C1
        sch['_R_C1'] = NULL_R
    else:
        sch['_R_C1'] = HUGE_R
        var_list += ['C1']

    if swcode & 16:  # D1
        sch['_R_D1'] = NULL_R
    else:
        sch['_R_D1'] = HUGE_R

    if swcode & 32:  # C3
        sch['_R_C3'] = NULL_R
    else:
        sch['_R_C3'] = HUGE_R
        var_list += ['C3']

    if swcode & 64:  # D3
        sch['_R_D3'] = NULL_R
    else:
        sch['_R_D3'] = HUGE_R

    if swcode & 128:  # C2
        sch['_R_C2'] = NULL_R
    else:
        sch['_R_C2'] = HUGE_R
        var_list += ['C2']

    session['Xi_variable'] = var_list


def Session_processAll(fileName='result.txt'):
    global FITTER_SUCCESS, VALUES_TOLERANCE, MAXFEV
    FITTER_SUCCESS = False
    ses_list = []
    best_ses = None
    best_misfit = 2  # заведомо большое число

    # создаем список сессий для старта
    for swcode in range(255):
        if not is_valid_switchers(swcode):
            continue

        code2 = 0
        next_code2 = True

        while next_code2:
            sch0 = Sch_init()
            ses = Session_create(sch0)
            next_code2 = Session_init_by_approximation(ses, swcode, code2, fileName)
            code2 += 1
            ses_list += [ses]

            if ses['misfit'] < IVCMP_TOLERANCE:  # условие останова удовлетворено
                best_ses = ses
                best_misfit = best_ses['misfit']
                print(ses['start_sch'])
                analysis_plot('FITTER SUCCESS')
                print('FITTER_SUCCESS!!\nmisfit = '+str(best_misfit))
                Sch_saveToFile(best_ses, fileName)
                print('good case!!')
                return

    # end_for
    print('pre init completed')
    # сортируем сессии, чтобы начать подбор с наиболее подходящих
    ses_list = sorted(ses_list, key=lambda s: s['misfit'])
    best_ses = ses_list[0]
    best_misfit = best_ses['misfit']

    # запускаем автоподбор, пока не будут удовлетворены условия останова
    for ses in ses_list:
        Session_run_fitter(ses)
        if (ses['misfit'] < best_misfit):
            best_misfit = ses['misfit']
            best_ses = ses
            print('misfit = '+str(best_misfit))
            if FITTER_SUCCESS:
                print(ses['result_sch'])
                analysis_plot('FITTER SUCCESS')

                print('FITTER_SUCCESS!!\nmisfit = '+str(best_misfit))
                Sch_saveToFile(best_ses, fileName)
                return
    # end_for

    # подбор завершился неудачно, выводим что есть
    print('FITTER routine unsuccessfull\nmisfit = '+str(best_ses['misfit']))
    Sch_saveToFile(best_ses, fileName)
    Session_run1(best_ses)
    analysis_plot('FITTER routine unsuccessfull')


def open_board(path):
    with open(path, "r") as dump_file:
        ivc_real = json.load(dump_file)
        return ivc_real
    return None

#############################################################################


def test2():
    sch = Sch_init()
    sch['R1'] = 1e2
    sch['C1'] = 1e-5
    sch['_R_C1'] = HUGE_R
    sch['R3'] = 1e3
    init_target_by_Sch(sch)
    print('test2()')
    Session_processAll('test2.txt')


def test3():
    sch = Sch_init()
    sch['R2'] = 1e2
    sch['_R_C2'] = HUGE_R
    sch['C2'] = 1e-7
    init_target_by_Sch(sch)
    print('test3()')
    Session_processAll('test3.txt')


def test4():
    sch = Sch_init()
    sch['R2'] = 1e2
    sch['_R_C2'] = HUGE_R
    sch['C2'] = 1e-7
    sch['R3'] = 1e3
    init_target_by_Sch(sch)
    print('test4()')
    Session_processAll('test4.txt')


def test5():
    sch = Sch_init()
    sch['R1'] = 1e2
    sch['C1'] = NONE_C
    sch['_R_C1'] = NULL_R
    sch['_R_D1'] = NULL_R

    sch['R2'] = NULL_R
    sch['_R_C2'] = HUGE_R
    sch['C2'] = 1e-7
    init_target_by_Sch(sch)
    print('test5()')
    Session_processAll('test5.txt')


def test_data_jsn(jsn_data, N, fileName='result.txt'):
    global Z123_sch
    gc.collect()
    print('\n')
    print(jsn_data)
    InitF, InitV, InitRcs, target_voltages, target_currents = init_target_from_jsnFile(jsn_data, N)
    init_target_Data(target_voltages, target_currents, initF=InitF, initV=InitV, initRcs=InitRcs, cycle=100)
    Session_processAll(fileName)


def test_circuit(circuitFile, resultFile='result.txt'):
    gc.collect()
    print('\n')
    print(circuitFile)
    init_target_by_circuitFile(circuitFile)
    Session_processAll(resultFile)
##############################################################


def main():
    k = 4
    test_data_jsn("vs_circuit_solver\\data\\100khz.json", k, 'vs_circuit_solver\\data\\100khz_{}.txt'.format(k))

if __name__ == '__main__':
    main()

    # Тестовый код с сырыми значениями
    # target_voltages = [-0.021364402829182096, -0.043362267221755214, -0.05730522931689329, -0.07996540177202412, -0.0986674560799491, -0.11657465624482724, -0.13006556774239905, -0.1473326646941493, -0.16391057305362058, -0.17962204107962565, -0.19736868696610627, -0.20735057697733458, -0.22039313161505417, -0.2343093011800344, -0.24943708368746734, -0.2569863007199858, -0.26621546117472616, -0.27481264061382005, -0.28169560720570935, -0.28905349632279775, -0.29006130686056375, -0.295348283232999, -0.29371794079640456, -0.3000439792012611, -0.29757401925836985, -0.3001002167318974, -0.29983801374372915, -0.2941089710229495, -0.2877997396765084, -0.28245344178733245, -0.2789915829431597, -0.2683187410245671, -0.26155951720218645, -0.2527694871631476, -0.23884815804125417, -0.2324553922050823, -0.21594136062993044, -0.20367704032562373, -0.18858820743250004, -0.17179172990075953, -0.1535157689829397, -0.13936412944678825, -0.12378160601076153, -0.10678805508593103, -0.084411584414049, -0.06950671174972531, -0.05172907544421448, -0.03245709328611624, -0.01563813189201769, 0.00855141249800883, 0.025377236954857364, 0.042214187055739476, 0.058029797017228685, 0.07837491348785751, 0.10086575073805157, 0.11462844092694037, 0.1286645217307816, 0.15012182404997979, 0.1651555306847408, 0.17829188549592417, 0.19762402954486374, 0.21190436303221005, 0.22295640873916342, 0.23224819202011157, 0.24515366921785567, 0.2543047596182349, 0.26506315022504173, 0.27284510467957257, 0.2802432364099473, 0.2853758558294593, 0.29164902193201175, 0.2975201317686522, 0.30068937938234735, 0.29999611079019667, 0.3036820406096977, 0.29878985990584056, 0.30115126767300154, 0.2927666997605613, 0.2876792626572526, 0.2812342159257204, 0.2808615326443627, 0.26645149835054127, 0.2609335095642022, 0.24981639262087169, 0.2419196210501858, 0.2276766690695221, 0.21694931122898764, 0.2013140644063413, 0.18718330426956803, 0.17096196735142744, 0.15727445671228474, 0.13848430861426547, 0.12162875877085182, 0.10439289402220209, 0.08815172124550069, 0.0685239269324982, 0.0498295530123106, 0.034503640936578525, 0.017240829840766872, 0.0012290560606675047]
    # target_currents = [-1.8454503097317802e-11, -1.781286341509498e-11, -1.7551925237930386e-11, -1.6945262037862598e-11, -1.7503501814285147e-11, -1.6500969343178112e-11, -1.6082969642810865e-11, -1.4700419035974535e-11, -1.4286202265675446e-11, -1.3869804608717821e-11, -1.302900175897418e-11, -1.1890819857272572e-11, -1.114571392946306e-11, -1.0617586356012384e-11, -1.0198561217790904e-11, -8.77004493981539e-12, -7.96604182164283e-12, -6.975003399575628e-12, -5.665584298164306e-12, -4.204257692684497e-12, -4.442123264530925e-12, -3.2195287221160576e-12, -2.904806569556428e-12, -4.657566618211492e-13, 5.506224550054051e-14, 2.665259654264664e-13, 1.718880020809096e-12, 2.6423844418856773e-12, 3.912541742090524e-12, 5.136608159918399e-12, 5.465641250822719e-12, 6.122866993221849e-12, 8.15792157593104e-12, 9.13334318174199e-12, 8.975053567035998e-12, 9.951333341985984e-12, 1.0532112685218612e-11, 1.172716356528274e-11, 1.3055610141776243e-11, 1.4586292793674017e-11, 1.4218437921991736e-11, 1.4752176464870328e-11, 1.6267279319310333e-11, 1.6871092725820307e-11, 1.6690413106810926e-11, 1.7732729581482633e-11, 1.827339955264207e-11, 1.8133918866977845e-11, 1.8714359360509475e-11, 1.8969909637359124e-11, 1.8924024478441412e-11, 1.928079366729875e-11, 1.9459822888899825e-11, 1.895435609498538e-11, 1.988851035451873e-11, 1.9124048612977558e-11, 1.9409611330652754e-11, 1.957058220346478e-11, 1.8791655020633328e-11, 1.874403917531808e-11, 1.9025854997325472e-11, 1.951034869703245e-11, 2.1484821889905793e-11, 2.3263615307753153e-11, 2.7165370828301867e-11, 3.280381211768725e-11, 3.991048750673726e-11, 4.927051349293826e-11, 6.102936394542807e-11, 7.340902805169316e-11, 8.642382997100632e-11, 9.696736369742505e-11, 1.0515189414233398e-10, 1.0939817083032466e-10, 1.0864410406419706e-10, 1.0381993552473133e-10, 9.378405665298594e-11, 8.087928981584292e-11, 6.623645972109099e-11, 5.082709299433179e-11, 3.587488871278757e-11, 2.3112531464458116e-11, 1.1967753883309884e-11, 3.0392632581121295e-12, -3.6360270995890154e-12, -8.794901052648776e-12, -1.1246921367657126e-11, -1.3694898225275867e-11, -1.5221693830522017e-11, -1.7447968574602116e-11, -1.790323786340546e-11, -1.778800691685112e-11, -1.8435299401831924e-11, -1.9107317421584717e-11, -1.9192326592264184e-11, -1.9654989529113503e-11, -1.9559453154220115e-11, -1.8738845401240422e-11, -1.871800656805242e-11, -1.9276334459417496e-11]
    #
    # # target_voltages = [-0.022303983693747302, -0.04258308847417828, -0.061552592687505306, -0.07984048032233755, -0.09674594739717157, -0.11381935393935792, -0.13336163514797839, -0.1496100939808471, -0.16575970334498394, -0.17808990296667193, -0.19537091823895222, -0.21055148781372976, -0.2241184371684157, -0.23663779710008784, -0.24583917463745913, -0.25484870698824, -0.2663842491982413, -0.2694305094939639, -0.278333875844606, -0.28573460917674365, -0.29100530182124124, -0.29417714906102005, -0.30057807390376695, -0.30082367834993035, -0.30042158119628337, -0.2970562814138667, -0.2962408519070912, -0.2942285968766555, -0.2874629296389486, -0.27973032769015477, -0.27783738561955684, -0.2720789093294497, -0.2586874964025822, -0.25041277632938796, -0.23918470859663266, -0.22590530893064784, -0.21224627731374937, -0.20208924005389944, -0.18387644095415334, -0.17350981324556053, -0.15441831341597353, -0.13786206897012046, -0.12194604775446659, -0.10557762035916124, -0.08748232260919744, -0.06945748556039072, -0.054158223601683875, -0.03173845233439013, -0.015284180520779697, 0.005178941303807875, 0.025864369565656183, 0.04217734659708809, 0.06032786990205169, 0.07923900587633263, 0.09760701308290372, 0.11214433892027152, 0.1333278592623079, 0.15006482103970376, 0.1648835947773955, 0.18262742080972016, 0.19789692560058697, 0.20844478601152167, 0.22231656549744813, 0.2321385885027514, 0.24931110352674152, 0.2548705307146573, 0.2675437277187994, 0.2789892068057761, 0.28106455463052954, 0.28534865001223514, 0.2923971773364995, 0.29794538939995674, 0.2961420286369324, 0.300988848418376, 0.2996923834274542, 0.2982322988081488, 0.2952928260584379, 0.2940868604208378, 0.2901016504929484, 0.2862743313087197, 0.27671355828042005, 0.2676199850716631, 0.26000142075948496, 0.25034552933353393, 0.24250448313281636, 0.22748893620747093, 0.2086584637259086, 0.19849646900846027, 0.19086887513632905, 0.1722464060322371, 0.15327001208832328, 0.14152908514963466, 0.12173565276495259, 0.10318423664120001, 0.08858736604880277, 0.07201022739917522, 0.05242994175527898, 0.03107639384970721, 0.01570945410536622, -0.0009162168777034703]
    # # target_currents = [-0.02577373440550276, -0.044634418525037974, -0.0641976316265999, -0.07624087575531002, -0.09741780128211774, -0.11938650279732778, -0.1363528543797237, -0.15143402414286028, -0.1658319461592371, -0.1801175425984708, -0.19491300922701527, -0.21060574341077237, -0.2233129029124872, -0.23549524022925342, -0.24725375265230265, -0.2538585185518085, -0.26372940813609086, -0.27391497283645005, -0.28080553090187843, -0.2885648830023235, -0.2927290592448105, -0.29153205367948426, -0.29926801969735994, -0.3044328342647475, -0.2976634281186189, -0.298389736533735, -0.2948429603591007, -0.29398729023949366, -0.28687410606545116, -0.28387570650855143, -0.2791084943099751, -0.26829888086281184, -0.2565691297489663, -0.25220239280077683, -0.23843542463303613, -0.22540546958560997, -0.21474891359907844, -0.19976849230461188, -0.1866048421893325, -0.17048607237131472, -0.1585789238146981, -0.14093072383405317, -0.12074124412645855, -0.10677535255391213, -0.08706242601662427, -0.07279449863883139, -0.05196986920328534, -0.028592994558459185, -0.012642697464379923, 0.004061944638470738, 0.022433103442623815, 0.04420849903033833, 0.06403039748089934, 0.08290896664139857, 0.09667697274817331, 0.11480535425096644, 0.1273679964674779, 0.14817787321342527, 0.16679743539834835, 0.18313342333977803, 0.19504151844886197, 0.2091843914275883, 0.2243069998859168, 0.23603232886378883, 0.24697057432682007, 0.25065679771790944, 0.26656050443302254, 0.27046094629883016, 0.27863906242062386, 0.2826804251918348, 0.29200869186259343, 0.29509068109010744, 0.30180783024997315, 0.30028190786206027, 0.2967265162926531, 0.3009333853704228, 0.29708625645264086, 0.29594412309963947, 0.2917620346353211, 0.28495450172071024, 0.2762880683744542, 0.2678276355205723, 0.2587215637122816, 0.25243149147953986, 0.23974728547792382, 0.228438444378606, 0.2176625264149819, 0.20222050058737057, 0.18243460117583205, 0.17182601836384248, 0.15653740802546054, 0.13947059017016708, 0.12136627899148945, 0.1024094345391997, 0.08389512818753565, 0.07070474032328146, 0.051088232887973865, 0.03296105003639654, 0.015495236268777648, 0.0037683169513853887]
    #
    #
    # InitF = 1000
    # InitV = 0.3
    # InitRcs = 0
    # fileName = 'spice4qucs'
    #
    # init_target_Data(target_voltages, target_currents, initF=InitF, initV=InitV, initRcs=InitRcs, cycle=100)
    # Session_processAll(fileName)

##############################################################################
