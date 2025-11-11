#!/bin/bash
################################################################################
# 批量提交冻结声子的VASP能带计算
#
# 用法：
#   bash batch_submit_vasp.sh <mode_number> <amplitudes>
#
# 示例：
#   bash batch_submit_vasp.sh 25 "0.000 0.050 0.100 0.150 0.200"
#
# 功能：
#   1. 为每个振幅创建独立目录
#   2. 复制INCAR, KPOINTS, POTCAR, POSCAR
#   3. 提交VASP任务（支持SLURM或直接运行）
#
################################################################################

# 检查参数
if [ $# -lt 2 ]; then
    echo "错误：参数不足"
    echo "用法：bash batch_submit_vasp.sh <mode_number> <amplitudes>"
    echo "示例：bash batch_submit_vasp.sh 25 \"0.000 0.050 0.100 0.150 0.200\""
    exit 1
fi

MODE_NUM=$1
AMPLITUDES=$2

# 配置文件路径（根据实际情况修改）
INCAR_TEMPLATE="../templates/INCAR_bands"
KPOINTS_TEMPLATE="../templates/KPOINTS_bands"
POTCAR_FILE="../POTCAR"
DISPLACED_DIR="./displaced_structures"

# VASP可执行文件
VASP_EXEC="vasp_std"  # 或 vasp_gam, vasp_ncl

# 并行设置
NPROCS=16  # MPI进程数

# 任务调度系统（slurm, pbs, 或 direct）
SCHEDULER="direct"  # 改为 "slurm" 或 "pbs" 如果使用调度器

################################################################################
# 主程序
################################################################################

echo "=========================================="
echo "批量提交冻结声子VASP计算"
echo "=========================================="
echo "模式编号：$MODE_NUM"
echo "振幅列表：$AMPLITUDES"
echo ""

# 检查模板文件
if [ ! -f "$INCAR_TEMPLATE" ]; then
    echo "错误：INCAR模板不存在 $INCAR_TEMPLATE"
    exit 1
fi

if [ ! -f "$KPOINTS_TEMPLATE" ]; then
    echo "错误：KPOINTS模板不存在 $KPOINTS_TEMPLATE"
    exit 1
fi

if [ ! -f "$POTCAR_FILE" ]; then
    echo "错误：POTCAR文件不存在 $POTCAR_FILE"
    exit 1
fi

# 对每个振幅创建目录并提交任务
for amp in $AMPLITUDES; do
    # 目录名
    DIR="mode${MODE_NUM}_amp${amp}"

    echo "----------------------------------------"
    echo "处理振幅：$amp Å"
    echo "目录：$DIR"

    # 创建目录
    mkdir -p "$DIR"
    cd "$DIR"

    # 复制输入文件
    echo "  复制输入文件..."
    cp "$INCAR_TEMPLATE" INCAR
    cp "$KPOINTS_TEMPLATE" KPOINTS
    cp "$POTCAR_FILE" POTCAR

    # 复制位移的POSCAR
    POSCAR_SRC="../${DISPLACED_DIR}/POSCAR_mode${MODE_NUM}_amp${amp}"
    if [ ! -f "$POSCAR_SRC" ]; then
        echo "  警告：POSCAR不存在 $POSCAR_SRC"
        echo "  跳过此振幅"
        cd ..
        continue
    fi
    cp "$POSCAR_SRC" POSCAR

    # 修改INCAR的SYSTEM标签（可选）
    sed -i "s/SYSTEM = .*/SYSTEM = Mode${MODE_NUM}_Amp${amp}/" INCAR

    # 提交任务
    echo "  提交VASP任务..."

    case $SCHEDULER in
        slurm)
            # SLURM提交脚本
            cat > submit.sh << EOF
#!/bin/bash
#SBATCH --job-name=mode${MODE_NUM}_${amp}
#SBATCH --nodes=1
#SBATCH --ntasks=${NPROCS}
#SBATCH --time=24:00:00
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err

module load vasp/5.4.4  # 根据实际环境修改

mpirun -np ${NPROCS} ${VASP_EXEC}
EOF
            chmod +x submit.sh
            sbatch submit.sh
            echo "  任务已提交（SLURM）: $(pwd)"
            ;;

        pbs)
            # PBS提交脚本
            cat > submit.sh << EOF
#!/bin/bash
#PBS -N mode${MODE_NUM}_${amp}
#PBS -l nodes=1:ppn=${NPROCS}
#PBS -l walltime=24:00:00
#PBS -o pbs.out
#PBS -e pbs.err

cd \$PBS_O_WORKDIR
module load vasp/5.4.4

mpirun -np ${NPROCS} ${VASP_EXEC}
EOF
            chmod +x submit.sh
            qsub submit.sh
            echo "  任务已提交（PBS）: $(pwd)"
            ;;

        direct)
            # 直接运行（适合小型计算或测试）
            echo "  直接运行VASP（后台）..."
            nohup mpirun -np ${NPROCS} ${VASP_EXEC} > vasp.log 2>&1 &
            PID=$!
            echo "  VASP进程ID: $PID"
            echo $PID > vasp.pid
            ;;

        *)
            echo "  错误：未知的调度器 $SCHEDULER"
            cd ..
            continue
            ;;
    esac

    # 返回上级目录
    cd ..

    # 短暂延迟（避免同时提交太多任务）
    sleep 1
done

echo ""
echo "=========================================="
echo "所有任务已提交完成"
echo "=========================================="
echo ""
echo "后续步骤："
echo "1. 等待所有计算完成"
echo "2. 检查每个目录的OUTCAR确认收敛"
echo "3. 运行分析脚本："
echo "   python analyze_band_splitting.py $MODE_NUM"
echo ""
