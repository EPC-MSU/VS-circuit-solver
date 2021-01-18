# vs_spice_solver
восстановление электронной схемы по ВАХ АЧХ, исследовательский проект



#### **Установка на Windows**

- Проверка работоспособности осуществлялась на стандартной версии питон python-3.6.8-amd64. Версия была установлена с атрибудами для всех пользователей и прописыванием путей в Path.
- Все зависимости прописаны в файле requirements.txt. Запуск установки зависимостей производится командой
  ***pip install -r requirements.txt***
- По мимо этих зависимостей необходимо скопировать на диск C:\Program Files\  папку  Spice64_dll. Если ее нет можно скачать https://sourceforge.net/projects/ngspice/files/ng-spice-rework/33/ файл ngspice-33_dll_64.zip.

#### **Запуск**

- Запустить этот проект (из корня): **python vs_circuit_solver\vs_circuit_solver.py** 
- Запуск осуществляется из директории расположения файла vs_circuit_solver.py командой **python vs_circuit_solver.py** .
- Файл **vs_circuit_solver.py** так же можно использовать как подключаемый модуль

#### Основные функции модуля и их использование

- **init_target_Data**(target_voltages, target_currents, initF=1e4, initV=5, initRcs=1e2, initSNR=120, cycle=10, ivcmpTolerance = 6e-2) - функция позволяет задать массивы тока и напряжения, а также параметры при которых происходило измерения если это необходимо.

- **Session_processAl**l(fileName='result.txt') - запуск процесса подбора схемы с сохранением в указанный файл.

- В текущем варианте в папке **VS_circuit_solver\vs_circuit_solver\data**   находятся  json файлы с данными для подбора схем. При независимом запуске приложения имя и номер записи в файле устанавливается в виде:

      k = 4
      test_data_jsn("data\\100khz.json",k,'data\\100khz_{}.txt'.format(k))
      
      или
      
      for k in range(10):
          test_data_jsn("data\\100khz.json",k,'data\\100khz_{}.txt'.format(k))

- Для запуска подбора по другим данным используется код в котором необходимо установить необходимые параметры.

      init_target_Data(target_voltages, target_currents, initF=InitF, initV=InitV, initRcs=InitRcs, cycle=100)
      Session_processAll(fileName)