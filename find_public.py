
import requests
from ecdsa import VerifyingKey, SECP256k1
import binascii


def get_public_key(address):
    # Get transaction history (replace with a different API if needed)
    url = f"https://blockchain.info/rawaddr/{address}"
    response = requests.get(url)
    data = response.json()

    if "txs" not in data or len(data["txs"]) == 0:
        print("No transactions found for this address.")
        return None

    # Find a transaction where the address was used as an INPUT
    for tx in data["txs"]:
        for inp in tx["inputs"]:
            if "prev_out" in inp and inp["prev_out"]["addr"] == address:
                txid = tx["hash"]  # Get the transaction ID
                print(f"Transaction found where address spent BTC: {txid}")
                
                # Fetch raw transaction details
                raw_tx_url = f"https://blockchain.info/rawtx/{txid}?format=json"
                raw_tx_response = requests.get(raw_tx_url)
                raw_tx_data = raw_tx_response.json()

                # Extract the public key from scriptSig
                for vin in raw_tx_data["inputs"]:
                    if "script" in vin and "prev_out" in vin and vin["prev_out"]["addr"] == address:
                        scriptSig = vin["script"]
                        print(f"Public Key Found: {scriptSig}")
                        return scriptSig

    print("No public key found. The address may not have spent any BTC.")
    return None

# Replace with an actual address
btc_address = "16jY7qLJnxb7CHZyqBP8qca9d51gAjyXQN"
public_key = get_public_key(btc_address)





if public_key:
    print("\nExtracted Public Key:", public_key)
    # Convert to compressed form
    print("Compressed Public Key:", public_key[-66:])
else:
    print("\nPublic key could not be retrieved.")
