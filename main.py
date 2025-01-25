# Добавление зависимостей для работы с FastAPI
from fastapi import Body, FastAPI, HTTPException
from fastapi.responses import FileResponse, JSONResponse
# Добавление зависимости для выполнения вычислений, связанных с диффурами
from scipy. integrate import odeint


# Инициализация приложения FastAPI
app = FastAPI()


# Вспомогательный класс F, выражающий полином третьей степени.
# Поля:
# a - коэффициент при x^3
# b - коэффициент при x^2
# c - коэффициент при x
# d - свободный коэффициент
# Методы:
# calc(x) - вычисляет значение полинома при заданном x
class F:
    def __init__(self, a, b, c, d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d


    def calc(self, x):
        return self.a* (x**3) + self.b * (x**2) + self.c * x + self.d


# Функция, выражающая систему уравнений Практической Работы 3
# Принимает:
# u - массив исследуемых параметров, представляющих собой функции mi(t), i=[1..33]
# t - массив временных точек от 0 до 1
# с - словарь нормировочных множителей вида {"Название": <Значение множителя>, ...}
# f - словарь полиномов и возмущений вида {"Название": <Экземпляр полинома>, ...}
# Возвращает:
# массив вычисленных выражений dmi/dt, i=[1..33]
def du3_dt(u,t,c,f):


        # Извлекаем элементы, которые являются исследуемыми mi(t), из массива u
        [m1_t, m2_t, m3_t, m4_t, m5_t, m6_t, m7_t, m8_t, m9_t, m10_t, m11_t, m12_t, m13_t, m14_t, m15_t, m16_t, m17_t,
        m18_t, m19_t, m20_t, m21_t, m22_t, m23_t, m24_t, m25_t, m26_t, m27_t, m28_t, m29_t, m30_t, m31_t, m32_t, m33_t] = u


        # Для каждого уравнения приводим запись формулы из задания. 
        # Для формулы:
        # dx1/dt = 1/normM1 * 
        # ( 
        #   f1(m3(t)) * f2(m4(t)) * f3(m15(t)) * f6(m19(t)) * f7(m22(t)) * (z1(t) + z2(t)) 
        #   -
        #   f8(m2(t)) * f9(m5(t)) * f10(m13(t)) * f11(m15(t)) * f12(m16(t)) * f13(m17(t)) * f14(m18(t)) *
        #   f15(m26(t)) * f16(m27(t)) * f17(m30(t)) * (z3(t))
        # )
        # программная запись примет вид:
        dm1_dt = (
            (1/c['normM1']) *
            (
                f['F1'].calc(m3_t) *
                f['F2'].calc(m4_t) *
                f['F3'].calc(m15_t) *
                f['F6'].calc(m19_t) *
                f['F7'].calc(m22_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z2'].calc(t)
                )
                -
                (
                    f['F8'].calc(m2_t) *
                    f['F9'].calc(m5_t) *
                    f['F10'].calc(m13_t) *
                    f['F11'].calc(m14_t) *
                    f['F12'].calc(m16_t) *
                    f['F13'].calc(m17_t) *
                    f['F14'].calc(m18_t) *
                    f['F15'].calc(m26_t) *
                    f['F16'].calc(m27_t) *
                    f['F17'].calc(m30_t)
                ) *
                (
                    f['Z3'].calc(t)
                )
            )
        )


        dm2_dt = (
            (1/c['normM2']) *
            (
                f['F18'].calc(m4_t) *
                f['F19'].calc(m8_t) *
                f['F20'].calc(m12_t) *
                f['F21'].calc(m13_t) *
                f['F22'].calc(m19_t) *
                f['F23'].calc(m26_t) *
                f['F24'].calc(m29_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z2'].calc(t) +
                    f['Z3'].calc(t) +
                    f['Z4'].calc(t) +
                    f['Z5'].calc(t) +
                    f['Z6'].calc(t) +
                    f['Z7'].calc(t)
                )
                -
                f['F30'].calc(m1_t) *
                f['F31'].calc(m3_t) *
                f['F32'].calc(m5_t)
            )
        )


        dm3_dt = (
            (1/c['normM3']) *
            (
                f['F33'].calc(m1_t) *
                f['F34'].calc(m2_t) *
                f['F35'].calc(m4_t) *
                f['F36'].calc(m5_t) *
                f['F37'].calc(m6_t) *
                f['F38'].calc(m8_t) *
                f['F39'].calc(m12_t) *
                f['F40'].calc(m15_t) *
                f['F41'].calc(m19_t) *
                (
                    f['Z2'].calc(t) +
                    f['Z6'].calc(t) +
                    f['Z7'].calc(t)
                )
                -
                f['F42'].calc(m13_t) *
                f['F43'].calc(m14_t) *
                f['F44'].calc(m16_t) *
                f['F45'].calc(m17_t) *
                f['F46'].calc(m22_t) *
                f['F47'].calc(m23_t) *
                f['F48'].calc(m24_t) *
                f['F49'].calc(m30_t)
            )
        )


        dm4_dt = (
            (1/c['normM4']) *
            (
                f['F50'].calc(m1_t) *
                f['F51'].calc(m2_t) *
                f['F52'].calc(m3_t) *
                f['F53'].calc(m7_t) *
                f['F54'].calc(m8_t) *
                f['F55'].calc(m12_t) *
                f['F56'].calc(m13_t) *
                f['F57'].calc(m17_t) *
                f['F58'].calc(m19_t) *
                f['F59'].calc(m22_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z2'].calc(t) +
                    f['Z6'].calc(t) +
                    f['Z7'].calc(t)
                )
                -
                f['F60'].calc(m23_t) *
                f['F61'].calc(m25_t) *
                f['F62'].calc(m32_t) *
                f['F63'].calc(m33_t)
            )
        )


        dm5_dt = (
            (1/c['normM5']) *
            (
                f['F64'].calc(m3_t) *
                f['F65'].calc(m7_t) *
                f['F66'].calc(m12_t) *
                f['F67'].calc(m15_t) *
                f['F68'].calc(m18_t) *
                f['F69'].calc(m27_t) *
                f['F70'].calc(m29_t) *
                (
                    f['Z2'].calc(t)
                )
                -
                f['F71'].calc(m22_t) *
                f['F72'].calc(m25_t) *
                f['F73'].calc(m30_t) *
                f['F74'].calc(m33_t) *
                (
                    f['Z3'].calc(t)
                )
            )
        )


        dm6_dt = (
            (1/c['normM6']) *
            (
                f['F75'].calc(m1_t) *
                f['F76'].calc(m3_t) *
                f['F77'].calc(m4_t) *
                f['F78'].calc(m12_t) *
                f['F79'].calc(m15_t) *
                f['F80'].calc(m18_t) *
                f['F81'].calc(m19_t) *
                f['F82'].calc(m33_t) *
                (
                    f['Z3'].calc(t) +
                    f['Z6'].calc(t)
                )
                -
                f['F83'].calc(m10_t) *
                f['F84'].calc(m13_t) *
                f['F85'].calc(m26_t) *
                f['F86'].calc(m30_t)
            )
        )


        dm7_dt = (
            (1/c['normM7']) *
            (
                f['F87'].calc(m1_t) *
                f['F88'].calc(m2_t) *
                f['F89'].calc(m8_t) *
                f['F90'].calc(m12_t) *
                f['F91'].calc(m13_t) *
                f['F92'].calc(m20_t) *
                f['F93'].calc(m21_t) *
                f['F94'].calc(m26_t) *
                f['F95'].calc(m30_t)
            ) *
                (
                    f['Z1'].calc(t) +
                    f['Z2'].calc(t) + #Z2
                    f['Z6'].calc(t) +
                    f['Z7'].calc(t)
                )
                -
                f['Z3'].calc(t) -
                f['Z4'].calc(t) -
                f['Z5'].calc(t)
        )


        dm8_dt = (
            (1/c['normM8']) *
            (
                f['F96'].calc(m2_t) *
                f['F97'].calc(m4_t) *
                f['F98'].calc(m7_t) *
                f['F99'].calc(m10_t) *
                f['F100'].calc(m12_t) *
                f['F101'].calc(m13_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z6'].calc(t)
                )
                -
                f['Z3'].calc(t) *
                f['F102'].calc(m20_t) *
                f['F103'].calc(m25_t) *
                f['F104'].calc(m28_t) *
                f['F105'].calc(m32_t) *
                f['F106'].calc(m33_t)
            )
        )


        dm9_dt = (
            (1/c['normM9']) *
            (
                f['F107'].calc(m1_t) *
                f['F108'].calc(m4_t) *
                f['F109'].calc(m7_t) *
                f['F110'].calc(m12_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z3'].calc(t)
                )
                -
                f['F111'].calc(m20_t) *
                f['F112'].calc(m28_t)
            )
        )


        dm10_dt = (
            (1/c['normM10']) *
            (
                f['F113'].calc(m3_t) *
                f['F114'].calc(m7_t) *
                f['F115'].calc(m12_t) *
                f['F116'].calc(m13_t) *
                f['F117'].calc(m18_t) *
                f['F119'].calc(m26_t) *
                (
                    f['Z3'].calc(t)
                )
                -
                f['F118'].calc(m20_t) *
                f['F120'].calc(m28_t)
            )
        )


        dm11_dt = (
            (1/c['normM11']) *
            (
                f['F121'].calc(m1_t) *
                f['F122'].calc(m2_t) *
                f['F123'].calc(m7_t) *
                f['F124'].calc(m8_t) *
                f['F125'].calc(m13_t) *
                f['F126'].calc(m23_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z3'].calc(t)
                )
                -
                f['F127'].calc(m20_t) *
                f['F128'].calc(m28_t)
            )
        )


        dm12_dt = (
            (1/c['normM12']) *
            (
                f['F129'].calc(m1_t) *
                f['F130'].calc(m2_t) *
                f['F131'].calc(m3_t) *
                f['F132'].calc(m4_t) *
                f['F133'].calc(m7_t) *
                f['F134'].calc(m8_t) *
                f['F135'].calc(m10_t) *
                f['F136'].calc(m26_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z3'].calc(t)
                )
                -
                f['F137'].calc(m28_t)
            )
        )


        dm13_dt = (
            (1/c['normM13']) *
            (
                f['F138'].calc(m2_t) *
                f['F139'].calc(m4_t) *
                f['F140'].calc(m7_t) *
                f['F141'].calc(m8_t) *
                f['F142'].calc(m10_t) *
                f['F143'].calc(m17_t) *
                f['F144'].calc(m18_t) *
                f['F145'].calc(m19_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z2'].calc(t) +
                    f['Z5'].calc(t) +
                    f['Z6'].calc(t)
                )
                -
                f['F146'].calc(m20_t) *
                f['F147'].calc(m28_t) *
                f['Z3'].calc(t)
            )
        )


        dm14_dt = (
            (1/c['normM14']) *
            (
                f['F281'].calc(m2_t) *
                f['F148'].calc(m3_t) *
                f['F149'].calc(m4_t) *
                f['F150'].calc(m7_t) *
                f['F151'].calc(m19_t) *
                f['F152'].calc(m26_t) *
                f['F153'].calc(m10_t) *
                f['F154'].calc(m27_t) *
                (
                    f['Z1'].calc(t)
                )
                -
                f['Z4'].calc(t)
            )
        )


        dm15_dt = (
            (1/c['normM15']) *
            (
                f['F155'].calc(m3_t) *
                f['F156'].calc(m4_t) *
                f['F157'].calc(m5_t) *
                f['F158'].calc(m6_t) *
                f['F159'].calc(m7_t) *
                f['F160'].calc(m16_t) *
                f['F161'].calc(m17_t) *
                f['F162'].calc(m18_t) *
                f['F163'].calc(m19_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z4'].calc(t)
                )
                -
                f['F164'].calc(m20_t) *
                f['F165'].calc(m33_t)
            )
        )


        dm16_dt = (
            (1/c['normM16']) *
            (
                f['F166'].calc(m3_t) *
                f['F167'].calc(m5_t) *
                f['F168'].calc(m7_t) *
                f['F169'].calc(m10_t) *
                f['F170'].calc(m12_t) *
                f['F171'].calc(m17_t) *
                f['F172'].calc(m18_t) *
                f['F173'].calc(m19_t) *
                f['F174'].calc(m24_t) *
                f['F175'].calc(m27_t) *
                f['F176'].calc(m30_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z6'].calc(t)
                )
                -
                f['F177'].calc(m20_t) *
                f['F178'].calc(m28_t) *
                f['F179'].calc(m33_t) *
                f['Z3'].calc(t)
            )
        )


        dm17_dt = (
            (1/c['normM17']) *
            (
                f['F283'].calc(m3_t) *
                f['F284'].calc(m4_t) *
                f['F285'].calc(m5_t) *
                f['F286'].calc(m16_t) *
                f['F287'].calc(m18_t) *
                f['F288'].calc(m19_t) *
                f['F289'].calc(m21_t) *
                f['F290'].calc(m27_t) *
                f['F291'].calc(m30_t) *
                (
                    f['Z2'].calc(t) +
                    f['Z7'].calc(t)
                )
            )
                -
                f['F292'].calc(m20_t) *
                f['F293'].calc(m28_t) *
                f['Z3'].calc(t)


        )


        dm18_dt = (
            (1/c['normM18']) *
            (
                f['F195'].calc(m1_t) *
                f['F196'].calc(m3_t) *
                f['F197'].calc(m4_t) *
                f['F198'].calc(m8_t) *
                f['F199'].calc(m9_t) *
                f['F200'].calc(m11_t) *
                f['F201'].calc(m16_t) *
                f['F202'].calc(m17_t) *
                f['F203'].calc(m19_t) *
                f['F204'].calc(m26_t) *
                f['F205'].calc(m30_t) *
                (
                    f['Z1'].calc(t)
                )
                -
                f['F206'].calc(m20_t) *
                f['F207'].calc(m28_t) *
                f['Z3'].calc(t)
            )
        )


        dm19_dt = (
            (1/c['normM19']) *
            (
                f['F208'].calc(m16_t) *
                f['F209'].calc(m5_t) *
                f['F210'].calc(m19_t) *
                (
                    f['Z2'].calc(t)
                )
                -
                f['Z3'].calc(t)
            )
        )


        dm20_dt = (
            (1/c['normM20']) *
            (
                f['Z3'].calc(t) -
                f['F211'].calc(m2_t) *
                f['F212'].calc(m4_t) *
                f['F213'].calc(m7_t) *
                f['F214'].calc(m8_t) *
                f['F215'].calc(m11_t) *
                f['F216'].calc(m13_t) *
                f['F217'].calc(m23_t) *
                f['F218'].calc(m25_t) *
                f['F219'].calc(m28_t) *
                f['F220'].calc(m31_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z6'].calc(t)
                )
            )
        )


        dm21_dt = (
            (1/c['normM21']) *
            (
                -f['F221'].calc(m4_t) *
                f['F222'].calc(m7_t) *
                f['F223'].calc(m8_t)
                +
                f['Z1'].calc(t) +
                f['Z6'].calc(t)
            )
        )


        dm22_dt = (
            (1/c['normM22']) *
            (
                f['Z1'].calc(t) +
                f['Z6'].calc(t) -
                f['F224'].calc(m6_t) *
                f['F225'].calc(m24_t)
            )
        )


        dm23_dt = (
            (1/c['normM23']) *
            (
                f['F226'].calc(m9_t) *
                f['F227'].calc(m11_t) *
                f['F228'].calc(m27_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z6'].calc(t)
                )
                -
                f['Z3'].calc(t)
            )
        )


        dm24_dt = (
            (1/c['normM24']) *
            (
                f['F229'].calc(m8_t) *
                f['F230'].calc(m11_t) *
                f['F231'].calc(m27_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z6'].calc(t)
                )
                -
                f['F232'].calc(m20_t) *
                f['F233'].calc(m28_t)
            )
        )


        dm25_dt = (
            (1/c['normM25']) *
            (
                f['F234'].calc(m9_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z6'].calc(t)
                )
                -
                f['F235'].calc(m26_t)
            )
        )


        dm26_dt = (
            (1/c['normM26']) *
            (
                (
                    f['Z1'].calc(t) +
                    f['Z6'].calc(t)
                )
                -
                f['F236'].calc(m13_t) *
                f['F237'].calc(m25_t)
            )
        )


        dm27_dt = (
            (1/c['normM27']) *
            (
                f['F238'].calc(m13_t) *
                f['F239'].calc(m14_t) *
                f['F240'].calc(m28_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z6'].calc(t)
                )
                -
                f['Z3'].calc(t)
            )
        )


        dm28_dt = (
            (1/c['normM28']) *
            (
                f['F241'].calc(m20_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z6'].calc(t)
                )
                -
                f['F242'].calc(m2_t) *
                f['F243'].calc(m4_t) *
                f['F245'].calc(m7_t) *
                f['F246'].calc(m8_t) *
                f['F247'].calc(m13_t) *
                f['F248'].calc(m19_t)
            )
        )


        dm29_dt = (
            (1/c['normM29']) *
            (
                f['F249'].calc(m13_t) *
                f['F250'].calc(m28_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z6'].calc(t)
                )
            )
        )


        dm30_dt = (
            (1/c['normM30']) *
            (
                f['F251'].calc(m17_t) *
                f['F252'].calc(m18_t) *
                (
                    f['Z2'].calc(t) +
                    f['Z3'].calc(t) +
                    f['Z7'].calc(t)
                )
            )
        )


        dm31_dt = (
            (1/c['normM31']) *
            (
                f['F253'].calc(m9_t) *
                f['F254'].calc(m32_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z6'].calc(t)
                )
            )
        )


        dm32_dt = (
            (1/c['normM32']) *
            (
                f['F255'].calc(m2_t) *
                f['F256'].calc(m4_t) *
                f['F257'].calc(m7_t) *
                f['F258'].calc(m11_t) *
                f['F259'].calc(m13_t) *
                f['F260'].calc(m14_t) *
                f['F261'].calc(m29_t) 
                -
                f['F262'].calc(m18_t) *
                f['F263'].calc(m20_t)
            )
        )


        dm33_dt = (
            (1/c['normM33']) *
            (
                f['F264'].calc(m20_t) *
                f['F265'].calc(m28_t) *
                (
                    f['Z1'].calc(t) +
                    f['Z2'].calc(t) + #Z2
                    f['Z6'].calc(t)
                ) 
                -
                f['F266'].calc(m5_t) *
                f['F267'].calc(m7_t) *
                f['F268'].calc(m8_t) *
                f['F269'].calc(m9_t) *
                f['F270'].calc(m10_t) *
                f['F271'].calc(m13_t) *
                f['F272'].calc(m14_t) *
                f['F273'].calc(m16_t) *
                f['F274'].calc(m17_t) *
                f['F275'].calc(m18_t) *
                f['F276'].calc(m19_t) *
                f['F277'].calc(m21_t) *
                f['F278'].calc(m23_t) *
                f['F279'].calc(m27_t)
            )
        )
        


        #Далее по аналогии со всеми уравнениями системы...


        # Возвращаем массив вычисленных dmi/dt, i=[1..33]
        return [dm1_dt, dm2_dt, dm3_dt, dm4_dt, dm5_dt, dm6_dt, dm7_dt, dm8_dt, dm9_dt, dm10_dt, dm11_dt, dm12_dt, dm13_dt, dm14_dt, dm15_dt, dm16_dt,
        dm17_dt, dm18_dt, dm19_dt, dm20_dt, dm21_dt, dm22_dt, dm23_dt, dm24_dt, dm25_dt, dm26_dt, dm27_dt, dm28_dt, dm29_dt, dm30_dt, dm31_dt, dm32_dt, dm33_dt]


# Указывает маршрут /lab3 до web-интерфеса, созданного для ввода значений из Практической Работы 3
@app.get("/lab3")
async def main3():
    # Возващаем в качетсве ответа HTML-страницу, расположенную по адресу /calculator_fastapi/lab3.html
    return FileResponse("./calculator_fastapi/lab3.html")


# Указываем маршрут /api/lab3/count, по которому можно получить результаты вычислений для Практической работы 3
# Тело данного POST-запрос имеет вид {"m0": <Словарь начальный значений>, "c": <Словарь нормировочных множителей>, "f": <Словарь полиномов и возмущений>}
# Записи <Словарь начальных значений> имеют вид "Имя значения": <значение>
# Записи <Словарь нормировочных множителей> имеют вид "Имя нормировочного множителя": <значение>
# Записи <Словарь полиномов и возмущений> имеют вид "Имя полинома (возмущения)": {"a": <значение>, "b": <значение>, "c": <значение>, "d": <значение> }
@app.post("/api/lab3/count")
async def count3(data  = Body()):
    try:
        # Извлекаем пришедшие начальные значения в t0
        t0 = data['m0']


        # Устанавливаем значения t в диапазоне от 0 до 1 с шагом 0.05, при которых будут проводиться вычисления
        t_span = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1]


        # Извлекаем пришедшие нормировочные множители в с
        c = data['c']


        # Извлекаем пришедшие коэффициенты полиномов и возмущений в uf
        uf = data['f']


        # На основе коэфициентов uf создадим удобные для работы экземпляры класса F, каждый из которых сопоставив с названием
        # Соберм пары "Имя полинома (возмущения)" : <Экземпляр класса F> в словарь f
        f = {}
        for key in uf:
            f[key] = F(uf[key]['a'],uf[key]['b'],uf[key]['c'],uf[key]['d'])


        # Передаем в качетсве аргументов в функцию для решения системы диффуров odeint 
        # Параметры:
        # du_dt - ранее описанная функция, выражающая систему уравнений Практической Работы 1
        # list(t0.values()) - список начальных значений
        # t_span - массив временных точек от 0 до 1
        # args=(c,f,) - дополнительный аргументы в виде констант и полиномов
        # full_output=True - для определения, решена ли система или нет
        # Получаем на выходе:
        # usolution - Решение системы в виде массива, элементы которого являются массивами решений для указанной точке на временной прямой
        # d - лог вычислений
        usolution,d = odeint(du3_dt, list(t0.values()), t_span, args=(c,f,), full_output=True)


        # Если решение не определено для заданных параметров - вернуть ошибку
        if (d["message"] != "Integration successful."):
            raise HTTPException(status_code=502)


        # Преобразуем значения для отправки их в качетсве JSON-ответа
        solution = [None] * len(usolution)
        idx = 0
        for elem in usolution:
            solution[idx] = list(elem)
            idx = idx + 1


        return JSONResponse(content={"message": solution})
    except Exception as e:
        print(f"An error occurred: {e}")
        # При некорректно предоставленных в теле запроса данных
        raise HTTPException(status_code=500)

