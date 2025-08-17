#!/usr/bin/env python3
"""
MDP (Membrane-Disrupting Peptides) 任务示例运行脚本
演示如何使用multiscale2运行CALVADOS MDP模拟
"""

import os
import sys
import subprocess

# 添加multiscale2模块到Python路径
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(current_dir))
multiscale2_path = os.path.join(project_root, 'multiscale2')
sys.path.insert(0, multiscale2_path)
print(f"添加multiscale2路径: {multiscale2_path}")
print(f"当前目录: {current_dir}")
print(f"项目根目录: {project_root}")

def check_conda_environment():
    """检查conda环境"""
    print("="*60)
    print("环境检查")
    print("="*60)
    
    # 检查是否在conda环境中
    conda_env = os.environ.get('CONDA_DEFAULT_ENV', 'base')
    print(f"当前conda环境: {conda_env}")
    
    if conda_env == 'base':
        print("⚠️  警告: 当前在base环境中")
        print("建议激活包含CALVADOS的conda环境")
    
    # 检查CALVADOS是否可用
    try:
        import calvados
        print("✓ CALVADOS模块可用")
    except ImportError:
        print("❌ CALVADOS模块不可用")
        print("请确保在包含CALVADOS的conda环境中运行")
        return False
    
    return True

def run_mdp_simulation():
    """运行MDP任务模拟"""
    print("="*60)
    print("MDP (Membrane-Disrupting Peptides) 任务示例")
    print("="*60)
    
    # 检查环境
    if not check_conda_environment():
        return
    
    # 设置环境
    os.environ['CUDA_VISIBLE_DEVICES'] = '0'
    print("✓ 设置CUDA_VISIBLE_DEVICES=0")
    
    # 导入生成器
    from calvados_generator import CalvadosGenerator
    
    # 初始化生成器
    config_path = "config.yaml"
    generator = CalvadosGenerator(config_path)
    
    print(f"✓ 生成器初始化成功")
    print(f"  任务类型: {generator.calvados_config.get('task_type', 'MDP')}")
    print(f"  模拟步数: {generator.calvados_config['steps']:,} 步")
    print(f"  输出频率: 每 {generator.calvados_config['wfreq']} 步")
    print(f"  温度: {generator.calvados_config['temp']} K")
    print(f"  离子浓度: {generator.calvados_config['ionic']} M")
    
    # 获取蛋白质配置信息
    protein_name = generator.protein_name
    protein_nmol = generator.protein_nmol
    print(f"  蛋白质名称: {protein_name}")
    print(f"  蛋白质数量: {protein_nmol}")
    
    # 设置输出目录
    output_dir = "output_mdp"
    
    print(f"\n准备生成CALVADOS脚本...")
    print(f"  输出目录: {output_dir}")
    print(f"  将自动生成components.yaml文件")
    
    # 生成并运行
    cmd = generator.generate_and_run(output_dir, protein_name, gpu_id=0, replica=1)
    
    print(f"✓ 脚本生成完成")
    print(f"  执行命令: {cmd}")
    
    # 询问是否运行
    response = input("\n是否运行CALVADOS模拟？(y/n): ")
    if response.lower() == 'y':
        print("\n开始运行CALVADOS模拟...")
        print("这可能需要几分钟时间...")
        
        # 设置环境变量
        env = os.environ.copy()
        env['PYTHONPATH'] = multiscale2_path + ':' + env.get('PYTHONPATH', '')
        
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, env=env)
        
        if result.returncode == 0:
            print("✓ CALVADOS prepare.py成功完成！")
            
            # 检查生成的CALVADOS目录
            calvados_dir = os.path.join(output_dir, f"{protein_name}_1")
            if os.path.exists(calvados_dir):
                print(f"\nCALVADOS输出目录: {calvados_dir}")
                print("生成的文件:")
                for file in os.listdir(calvados_dir):
                    file_path = os.path.join(calvados_dir, file)
                    if os.path.isfile(file_path):
                        size = os.path.getsize(file_path)
                        print(f"  {file}: {size:,} bytes")
                
                # 运行实际的CALVADOS模拟
                run_script = os.path.join(calvados_dir, 'run.py')
                if os.path.exists(run_script):
                    print(f"\n运行CALVADOS模拟: {run_script}")
                    run_result = subprocess.run(f"cd {calvados_dir} && python run.py", 
                                              shell=True, capture_output=True, text=True)
                    
                    if run_result.returncode == 0:
                        print("✓ CALVADOS模拟成功完成！")
                        
                        # 检查DCD轨迹文件
                        dcd_file = os.path.join(calvados_dir, f"{protein_name}_1.dcd")
                        if os.path.exists(dcd_file):
                            size = os.path.getsize(dcd_file)
                            print(f"\n✓ 找到DCD轨迹文件: {dcd_file}")
                            print(f"  文件大小: {size:,} bytes")
                            
                            if size > 100000:  # 大于100KB
                                print("✓ DCD文件包含轨迹数据")
                            else:
                                print("⚠️  DCD文件可能较小")
                        else:
                            print(f"\n❌ 未找到DCD轨迹文件: {dcd_file}")
                        
                        # 列出所有最终文件
                        print(f"\n所有最终文件:")
                        for file in os.listdir(calvados_dir):
                            file_path = os.path.join(calvados_dir, file)
                            if os.path.isfile(file_path):
                                size = os.path.getsize(file_path)
                                print(f"  {file}: {size:,} bytes")
                        
                    else:
                        print("❌ CALVADOS模拟失败")
                        print("\n错误输出:")
                        print(run_result.stderr)
                        print("\n标准输出:")
                        print(run_result.stdout)
                else:
                    print("❌ 未找到run.py文件")
            else:
                print("❌ 未找到CALVADOS输出目录")
        else:
            print("❌ CALVADOS prepare.py失败")
            print("\n错误输出:")
            print(result.stderr)
            print("\n标准输出:")
            print(result.stdout)
    else:
        print("跳过运行，仅生成脚本文件")
        print(f"你可以手动运行: {cmd}")

def main():
    """主函数"""
    print("MDP任务示例")
    print("="*60)
    print("这个示例演示如何运行MDP (Membrane-Disrupting Peptides) 任务")
    print("使用CALVADOS进行粗粒化分子动力学模拟")
    print()
    print("⚠️  重要提示:")
    print("请确保已激活包含CALVADOS的conda环境")
    print("例如: conda activate your_calvados_env")
    print()
    
    run_mdp_simulation()
    
    print("\n" + "="*60)
    print("MDP任务示例完成！")
    print("="*60)

if __name__ == "__main__":
    main()
