#!/usr/bin/env python3
"""
测试MultiscaleWorkflow的配置读取功能
"""

import sys
import os
import logging

# 添加项目根目录到Python路径
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from multiscale2.workflow import MultiscaleWorkflow, ConfigurationError

def test_config_loading():
    """测试配置加载功能"""
    print("="*60)
    print("测试配置加载功能")
    print("="*60)
    
    try:
        # 配置文件路径
        config_path = "00_input/config.yaml"
        
        # 检查配置文件是否存在
        if not os.path.exists(config_path):
            print(f"❌ 配置文件不存在: {config_path}")
            return False
        
        print(f"✅ 配置文件存在: {config_path}")
        
        # 创建MultiscaleWorkflow实例
        workflow = MultiscaleWorkflow(config_path)
        
        print(f"✅ 成功创建MultiscaleWorkflow实例")
        print(f"   项目根目录: {workflow.project_root}")
        print(f"   项目名称: {workflow.config.get('project_name', 'Unknown')}")
        
        # 显示配置信息
        print("\n📋 配置信息摘要:")
        print(f"   - 项目名称: {workflow.config.get('project_name')}")
        print(f"   - 温度: {workflow.config.get('cg_calvados', {}).get('temp')} K")
        print(f"   - 模拟步数: {workflow.config.get('cg_calvados', {}).get('steps')}")
        print(f"   - 盒子尺寸: {workflow.config.get('cg_calvados', {}).get('box')}")
        print(f"   - 平台: {workflow.config.get('cg_calvados', {}).get('platform')}")
        
        return True
        
    except ConfigurationError as e:
        print(f"❌ 配置错误: {e}")
        return False
    except Exception as e:
        print(f"❌ 未知错误: {e}")
        return False

def test_input_files():
    """测试输入文件验证"""
    print("\n" + "="*60)
    print("测试输入文件验证")
    print("="*60)
    
    try:
        config_path = "00_input/config.yaml"
        workflow = MultiscaleWorkflow(config_path)
        
        # 检查输入文件
        input_files = workflow.config.get('input_files', {})
        
        for file_key, file_path in input_files.items():
            full_path = workflow._get_path(file_path)
            if os.path.exists(full_path):
                print(f"✅ {file_key}: {full_path}")
            else:
                print(f"❌ {file_key}: {full_path} (文件不存在)")
                return False
        
        return True
        
    except Exception as e:
        print(f"❌ 输入文件验证失败: {e}")
        return False

def test_components_file():
    """测试components文件验证"""
    print("\n" + "="*60)
    print("测试components文件验证")
    print("="*60)
    
    try:
        config_path = "00_input/config.yaml"
        workflow = MultiscaleWorkflow(config_path)
        
        # 验证components文件
        components_path = workflow._validate_components_file()
        print(f"✅ Components文件验证成功: {components_path}")
        
        # 读取并显示components信息
        import yaml
        with open(components_path, 'r') as f:
            components = yaml.safe_load(f)
        
        protein_config = components.get('system', {}).get('protein', {})
        print(f"   - 蛋白质分子数: {protein_config.get('nmol')}")
        print(f"   - FASTA文件: {protein_config.get('ffasta')}")
        
        return True
        
    except Exception as e:
        print(f"❌ Components文件验证失败: {e}")
        return False

def test_calvados_config():
    """测试CALVADOS配置验证"""
    print("\n" + "="*60)
    print("测试CALVADOS配置验证")
    print("="*60)
    
    try:
        config_path = "00_input/config.yaml"
        workflow = MultiscaleWorkflow(config_path)
        
        # 获取CALVADOS配置
        cg_config = workflow.config.get('cg_calvados', {})
        
        print("📋 CALVADOS配置参数:")
        for key, value in cg_config.items():
            print(f"   - {key}: {value}")
        
        print("\n✅ CALVADOS配置验证通过")
        return True
        
    except Exception as e:
        print(f"❌ CALVADOS配置验证失败: {e}")
        return False

def main():
    """主测试函数"""
    print("🚀 开始测试MultiscaleWorkflow配置读取功能")
    
    # 设置日志级别
    logging.basicConfig(level=logging.INFO)
    
    tests = [
        ("配置加载", test_config_loading),
        ("输入文件验证", test_input_files),
        ("Components文件验证", test_components_file),
        ("CALVADOS配置验证", test_calvados_config),
    ]
    
    passed = 0
    total = len(tests)
    
    for test_name, test_func in tests:
        print(f"\n🧪 运行测试: {test_name}")
        if test_func():
            passed += 1
            print(f"✅ {test_name} 测试通过")
        else:
            print(f"❌ {test_name} 测试失败")
    
    print("\n" + "="*60)
    print(f"测试结果: {passed}/{total} 通过")
    
    if passed == total:
        print("🎉 所有测试通过！配置读取功能正常工作。")
        return True
    else:
        print("⚠️  部分测试失败，请检查配置文件和输入文件。")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
