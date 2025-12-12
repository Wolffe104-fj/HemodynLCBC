# HemodynLCBC
Pet-project that simulates the hemodynamics of the large circle of blood circulation. 
The project was inspired by a real-life task for medical engineers that was presented on [the website of the Sechenov MSMU.](https://theranostic.sechenov.ru/caseheart)

## 📁 Content

<details>
<summary>Original task and the problems associated with it</summary>

  <br>
  
  You can find the original task [here](https://theranostic.sechenov.ru/caseheart). Here, I will translate its text into English and explain what difficulties arose while trying to solve the task.
  
  <img align="right" width="35%" height="500" alt="OriginalTask_EN" src="https://github.com/user-attachments/assets/c4312afc-cf27-4dcb-a73d-0f69d3bffbf6" />
  
  Models of the circulatory system are used to study the processes occurring within the circulatory system and the effects of medical devices such as artificial heart valves and auxiliary circulation devices. These models allow for personalized determination of the optimal patient condition and can be integrated into medical decision support systems.

It is necessary to simulate the large circle of blood circulation using software packages such as Matlab, Python, or C++. The system of differential equations for this simulation is provided in the appendix, along with the Euler method for numerical solution. The model parameters can be found in the appendix table.

In the appendix, you will find a brief overview of the theory and formulas required for modeling, as well as the model parameters and initial values for the system of differential equations.

The answer to the question is the value of the systolic pressure in the aorta during the last cycle, rounded to the nearest integer.

  <br>
  
  The main challenge I faced was a lack of materials necessary to solve the task. Therefore, I decided to approach it from a different angle. I reviewed the current solution to the issue and used AI to gather relevant initial data. This process helped me determine the approximate initial parameters for the model:

| Notation          | Value    | Unit of Measure   | Brief Description                               |
|       :---:       |   :---:  |  :---:            |  :---                                           |
| R1                | 1        | mmHg * s / mL     | Vascular segment resistance (Between C2 and C3) |
| R2                | 0.005    | mmHg * s / mL     | Mitral valve resistance                         |
| R3                | 0.013    | mmHg * s / mL     | Aortic valve resistance                         |
| R4                | 0.0398   | mmHg * s / mL     | Peripheral vascular resistance                  |
| C2                | 4.4      | mL / mmHg         | Capacitance of the first arterial chamber       |
| C3                | 1.33     | mL / mmHg         | Capacitance of the arterial chamber             |
| C4                | 0.8      | mL / mmHg         | Aortic capacitance                              |
| L                 | 0.0005   | mmHg * s**2 / mL  | Inductance (blood flow inertia)                 |
| dt                | 0.01     | s                 | Integration step (Euler's method)               |
| HR                | 75       | beats/min         | Heart rate                                      |
| Umax              | 2        | conv. units       | Maximum ventricular elastance                   |
| Umin⁡              | 0.05     | conv. units       | Minimum ventricular elastance                   |

And also, the X vector:

| Notation         | Value    | Unit of Measure | Brief Description     |
|  :---:           |  :---:   |  :---:          |  :---                 |
| x1               | 8        | mmHg            | Ventricular pressure  |
| x2               | 7.3      | mmHg            | Atrial pressure       |
| x3               | 70       | mmHg            | Arterial pressure     |
| x4               | 75       | mmHg            | Aortic pressure       |
| x5               | 20       | mL/s            | Blood flow velocity   |

As well as formulas. Compact recording of the ODE system for modeling a large circle of blood circulation:
<img width="60%" height="60%" alt="ODE_1" src="https://github.com/user-attachments/assets/85aafd8e-674e-4aa3-95bc-a23563cac3dd" />
<img width="60%" height="60%" alt="ODE_2" src="https://github.com/user-attachments/assets/b1ab7974-5163-4802-8372-b71afe70406a" />

Now you can use all of this to solve the task from the beginning. 
If anything, I have tried to comment on my code as much as possible, so that you can understand it in case of any difficulties.

### ⚠️ Important:
1. Some of the comments and designations in this draft may be incorrect, as it was developed and studied independently, and I am not an expert in this field.
2. This project cannot be used to address real-world medical issues, as it is a learning project and does not fully represent what happens in the human body.

I'm always open to communication and constructive criticism. If you have more expertise on the subject — please just write to me about the mistake, I'd be very grateful. As the saying goes, he who makes no mistakes, makes nothing.

<br>

> P.S. Here are some useful materials that helped me complete the task. I hope they will be helpful to you too:
>
> P. I. Begun - Biomechanics (ISBN 5-7325-0309-5)
> > There are chapters on the biomechanics of the heart and the vascular system. Even during my self-study of biomechanics, I came across this book. It is written well and in a simple language.
> 
> B. I. Tkachenko - The basis of human physiology. A manual for the higher educational schools, in 2 volumes (ISBN 5-86050-055-6)
> > There are also chapters on hemodynamics. For the most part, I've already used this tutorial to test my work. It also contains useful materials for completing tasks.
> 
> L. Formaggia, A. Quarteroni, A. Veneziani Eds. - Cardiovascular mathematics. Modeling and simulation of the circulatory system (ISBN 978-88-470-1151-9)
> > The book is very close to the topic, as it describes it in detail. However, it was difficult for me to read and understand. There are a lot of formulas involved. Nevertheless, the book is still useful and may be more helpful to you than it was to me.
<br>
</details>



<details>
<summary>My decision and the results ⚠️[Writing in progress]⚠️</summary>
  
  <img align="right" width="55%"  alt="Figure_1 1" src="https://github.com/user-attachments/assets/815ddf18-dcda-4fdf-b81c-f1a054d41da7" />
  
  <img align="left" width="55%" alt="Figure_2 1" src="https://github.com/user-attachments/assets/80fa45a4-75c6-41d3-b1fc-d6c42dc3bde8" />
  
  <img align="right" width="55%" alt="Figure_3 1" src="https://github.com/user-attachments/assets/4a9964cf-afd5-40b3-8fa9-d2e9d25988e2" />
  
  <img align="left" width="55%" alt="Figure_4 1" src="https://github.com/user-attachments/assets/50bcd109-db26-4fbb-8e20-9808f4648240" />

<br>
<br>

<br>
</details>

## 📌 Quick Start & Installation
### 📋 Requirements
All dependencies are pinned in requirements.txt:
| Package         | Version     |
| --------------- | ----------- |
| numpy           | 2.3.5       |
| matplotlib      | 3.10.7      |
| scipy           | 1.16.3      |
| contourpy       | 1.3.3       |
| cycler          | 0.12.1      |
| fonttools       | 4.61.0      |
| kiwisolver      | 1.4.9       |
| packaging       | 25.0        |
| pillow          | 12.0.0      |
| pyparsing       | 3.2.5       |
| python-dateutil | 2.9.0.post0 |
| six             | 1.17.0      |

### 🛠️ Installation
1. Clone repository (All OS)
```
git clone https://github.com/Wolffe104-fj/HemodynLCBC.git
cd HemodynLCBC
```

2. Create & activate virtual environment
```
python -m venv venv
```

| Windows               | Linux/macOS              |
| --------------------- | ------------------------ |
| venv\\Scripts\\activate | source venv/bin/activate |

3. Install dependencies (All OS)
```
pip install -r requirements.txt
```

4. Run simulation (All OS)
```
python main.py
```
