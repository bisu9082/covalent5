# covalent5 — 우분투 실행 가이드 (폴더·명령 정확판)

> 규칙: 각 단계마다 **[폴더]** 에서 **[명령]** 실행 → **[생성물]** 확인.
> 맨 처음 `COV5` 경로만 한 번 정하면, 이후 명령은 어느 위치에서 실행해도 됩니다.

---

## ★ STEP 0. 경로 설정 (터미널 열고 제일 먼저, 1회)

폴더가 어디 있느냐에 따라 한 줄만 고르세요.

```bash
# (A) Windows의 Downloads를 WSL 우분투에서 쓰는 경우:
export COV5="/mnt/c/Users/kkaan/Downloads/claude_research/covalent5"

# (B) 별도 리눅스 머신으로 폴더를 복사해 온 경우 (예: 홈에 복사):
# export COV5="$HOME/covalent5"

cd "$COV5" && pwd        # 출력이 covalent5 폴더 경로면 OK
ls                       # build_molecules.py, data, dft, docking, latex 가 보여야 함
```
이후 모든 단계는 이 터미널에서 진행(`COV5` 변수 유지).

---

## STEP 1. 환경 준비 (1회만)

**[폴더] 아무데나** (설치 작업)
```bash
conda create -n covalent5 python=3.10 -y
conda activate covalent5
pip install rdkit
conda install -c conda-forge openbabel -y
```
별도 수동 설치(스크립트로 불가): **ORCA 6.x**, **GNINA**.
확인:
```bash
which orca      # 경로 나오면 OK
gnina --version # 버전 나오면 OK
```

---

## STEP 2. 분자집합·입력 클린 재생성 (가벼움, 수 초)

**[폴더] `$COV5`**
```bash
cd "$COV5"
rm -rf data/structures dft docking/ligands docking/out_noncovalent docking/out_covalent
python build_molecules.py
python generate_dft_inputs.py --nprocs 8 --maxcore 3500     # nprocs=실제 코어수
python docking/prepare_ligands.py
```
**[생성물]**
- `data/structures/<code>/*.xyz` (30개) + `data/molecule_manifest.csv`
- `dft/<code>/<id>/{opt.inp,sp.inp,run.sh}` (30잡) + `dft/job_manifest.csv` + `dft/run_all.sh`
- `docking/ligands/{noncovalent,covalent}/*.sdf` + `docking/ligands/ligand_manifest.csv`

확인:
```bash
wc -l data/molecule_manifest.csv          # 31 (헤더+30)
find dft -name opt.inp | wc -l            # 30
cat docking/ligands/ligand_manifest.csv | wc -l   # 31
```

---

## STEP 3. DFT 테스트 1개 (먼저 반드시 — 환경 검증)

**[폴더] `$COV5/dft/DDVP/DDVP_00_achiral`**
```bash
cd "$COV5/dft/DDVP/DDVP_00_achiral"
orca opt.inp > opt.out 2>&1
orca sp.inp  > sp.out  2>&1
```
**[확인]**
```bash
ls opt.xyz                                 # opt 결과 좌표 생겼는지
grep "FINAL SINGLE POINT ENERGY" sp.out    # 에너지 한 줄 나오면 성공
grep "HURRAY"  opt.out                      # 최적화 수렴 메시지
```
실패하면 여기서 멈추고 에러를 나에게 보내주세요 (ORCA 설치/키워드 문제).

---

## STEP 4. DFT 전체 실행 (★무거움, 수 시간~)

**[폴더] `$COV5`**
```bash
cd "$COV5"
screen -S dft        # 또는 nohup. 세션 끊겨도 계속 돌게
bash dft/run_all.sh
# (screen에서 빠져나오기: Ctrl+A 누른 뒤 D)
```
**[생성물]** 각 잡 폴더에 `opt.out`, `sp.out`, `opt.xyz`.
진행 확인:
```bash
ls dft/*/*/sp.out | wc -l                  # 완료된 잡 수 (목표 30)
```

---

## STEP 5. 수용체 준비 (가벼움)

**[폴더] `$COV5/docking`**
```bash
cd "$COV5/docking"
bash prepare_receptor.sh
```
**[생성물]** `receptor_H.pdb`, `receptor.pdbqt`

---

## STEP 6. 도킹 실행 (GPU)

**[폴더] `$COV5/docking`**
```bash
cd "$COV5/docking"
bash run_gnina_noncovalent.sh     # 비공유 30 → out_noncovalent/
bash run_gnina_covalent.sh        # 공유 26 → out_covalent/
python parse_docking.py           # → docking_results.csv
```
**[확인]**
```bash
cat docking_results.csv | head     # affinity/CNNscore 표
```
**검증:** VX(6CQT/6CQX) 공유 포즈가 실험 부가물 기하 재현하는지 (docking/README.md).

---

## 끝나면 나에게 보낼 것
1. `dft/job_manifest.csv` + 각 `sp.out` (또는 압축) — 디스크립터 추출
2. `docking/docking_results.csv` — RQ1 분석
→ 받으면 디스크립터 추출 + RQ1(비공유 vs 공유 순위) 분석으로 진행.

---
### 한눈에
| STEP | 폴더 | 명령 | 비고 |
|---|---|---|---|
| 0 | (터미널) | `export COV5=...; cd $COV5` | 1회 |
| 1 | 아무데나 | conda/pip/ORCA/GNINA 설치 | 1회 |
| 2 | `$COV5` | `python build_molecules.py` 등 3개 | 가벼움 |
| 3 | `$COV5/dft/DDVP/DDVP_00_achiral` | `orca opt.inp...` | 테스트 |
| 4 | `$COV5` | `bash dft/run_all.sh` | ★무거움 |
| 5 | `$COV5/docking` | `bash prepare_receptor.sh` | 가벼움 |
| 6 | `$COV5/docking` | `bash run_gnina_*.sh` | GPU |
