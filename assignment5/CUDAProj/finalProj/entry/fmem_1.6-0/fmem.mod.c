#include <linux/module.h>
#include <linux/vermagic.h>
#include <linux/compiler.h>

MODULE_INFO(vermagic, VERMAGIC_STRING);

__visible struct module __this_module
__attribute__((section(".gnu.linkonce.this_module"))) = {
	.name = KBUILD_MODNAME,
	.init = init_module,
#ifdef CONFIG_MODULE_UNLOAD
	.exit = cleanup_module,
#endif
	.arch = MODULE_ARCH_INIT,
};

static const struct modversion_info ____versions[]
__used
__attribute__((section("__versions"))) = {
	{ 0x9412fa01, __VMLINUX_SYMBOL_STR(module_layout) },
	{ 0x6bc3fbc0, __VMLINUX_SYMBOL_STR(__unregister_chrdev) },
	{ 0xc9296d71, __VMLINUX_SYMBOL_STR(device_destroy) },
	{ 0x4a1dbd04, __VMLINUX_SYMBOL_STR(__register_chrdev) },
	{ 0x1eecf3fc, __VMLINUX_SYMBOL_STR(mutex_unlock) },
	{ 0x4f8b5ddb, __VMLINUX_SYMBOL_STR(_copy_to_user) },
	{ 0x27e1a049, __VMLINUX_SYMBOL_STR(printk) },
	{ 0x91e14f04, __VMLINUX_SYMBOL_STR(mutex_lock) },
	{ 0xda6982e5, __VMLINUX_SYMBOL_STR(device_create) },
	{ 0xc6cbbc89, __VMLINUX_SYMBOL_STR(capable) },
	{ 0x42c8de35, __VMLINUX_SYMBOL_STR(ioremap_nocache) },
	{ 0xbdfb6dbb, __VMLINUX_SYMBOL_STR(__fentry__) },
	{ 0xedc03953, __VMLINUX_SYMBOL_STR(iounmap) },
	{ 0x268f87ce, __VMLINUX_SYMBOL_STR(class_destroy) },
	{ 0xf6e19a00, __VMLINUX_SYMBOL_STR(directly_mappable_cdev_bdi) },
	{ 0x8d2e268b, __VMLINUX_SYMBOL_STR(param_ops_ulong) },
	{ 0xd7bd463f, __VMLINUX_SYMBOL_STR(__class_create) },
};

static const char __module_depends[]
__used
__attribute__((section(".modinfo"))) =
"depends=";


MODULE_INFO(srcversion, "3C96174E0257438395AE052");
